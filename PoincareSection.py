import numpy as np

class Plane:
    def __init__(self, normal, offset):
        self.offset = offset
        self.normal = normal / np.linalg.norm(normal)

    def WhichSide(self, A):
        assert  A.shape == self.normal.shape

        dot = np.dot(A, self.normal)
        if dot > self.offset:
            return 1
        elif dot < self.offset:
            return -1
        else:
            return 0

    def CheckIntersect(self, A, B):
        sideA = self.WhichSide(A)
        sideB = self.WhichSide(B)

        return sideA * sideB <= 0, sideA, sideB

    def GetIntersection(self, A, B):
        component = np.dot(B-A, self.normal)

        if np.abs(component) < 1e-8:
            return None
        dist = self.Distance(A)
        t = dist/component
        return A + t*(B-A)

    def Distance(self, A):
        return np.abs(np.dot(A, self.normal) - self.offset)

class PoincareSection:
    def __init__(self, plane, integrator=None, N=None, problem=None, t0=None, tF=None):
        self.plane = plane
        self.problem = problem
        self.t0 = t0
        self.tF = tF
        self.integrator = integrator
        self.N = N

    def ComputeFromIntegration(self, x0):
        intersections = []
        assert integrator is not None and problem is not None
        generator = self.integrator.IntegrateLazy(self.problem, self.t0, self.tF, x0, self.N)

        oldstate = next(generator)
        step = 1
        for state in generator:
            step += 1
            theyInters, direction, _ = self.plane.CheckIntersect(oldstate, state)
            if theyInters and direction == -1:
                inters = self.plane.GetIntersection(oldstate, state)
                intersections.append(inters)

            oldstate = state

        return intersections

    def ComputeFromOrbit(self, orbit):
        intersections = []
        
        oldstate = orbit[0,:]
        for i in range(1, orbit.shape[0]):
            state = orbit[i,:]
            theyInters, direction, _ = self.plane.CheckIntersect(oldstate, state)
            if theyInters and direction == -1:
                inters = self.plane.GetIntersection(oldstate, state)
                intersections.append(inters)

            oldstate = state

        return intersections

    def ComputeFromOrbitAxis(self, orbit):
        offset = self.plane.offset 
        assert np.all(self.plane.normal >=0) and np.sum(self.plane.normal == 1)
        coordinate = np.nonzero(self.plane.normal)[0]

        shift = np.zeros_like(self.plane.normal)
        shift[coordinate] = offset

        orbit_shifted = orbit - shift
        crossings_idx = np.where(np.logical_and(orbit_shifted[0:-1, coordinate] <= 0, orbit_shifted[1:, coordinate] >= 0))[0]

        intersections = [None] * len(crossings_idx)


        for i,idx in enumerate(crossings_idx):
            inters = self.plane.GetIntersection(orbit[idx,:], orbit[idx+1,:])
            intersections[i] = inters

        return intersections


    def ComputeFromOrbitAxisTimes(self, orbit, times):
        offset = self.plane.offset 
        assert np.all(self.plane.normal >=0) and np.sum(self.plane.normal == 1)
        coordinate = np.nonzero(self.plane.normal)[0]

        shift = np.zeros_like(self.plane.normal)
        shift[coordinate] = offset

        orbit_shifted = orbit - shift
        crossings_idx = np.where(np.logical_and(orbit_shifted[0:-1, coordinate] <= 0, orbit_shifted[1:, coordinate] >= 0))[0]
        crossings_times = times[crossings_idx]

        intersections = [None] * len(crossings_idx)


        for i,idx in enumerate(crossings_idx):
            inters = self.plane.GetIntersection(orbit[idx,:], orbit[idx+1,:])
            intersections[i] = inters

        return intersections, crossings_times

