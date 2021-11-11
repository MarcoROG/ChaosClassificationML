import numpy as np

class Corrector:
    def __init__(self, constant, eps_ord, h_ord):
        self.c = constant
        self.h_ord = h_ord
        self.eps_ord = eps_ord

class SplittingMethod:
    def __init__(self, cs, ds, corrector=None):
        self.cs = cs
        self.ds = ds
        assert len(cs) == len(ds)

        self.corrector = corrector

    def IntegrateLazy(self, problem, t0, tF, x0, N):
        FluxA = problem.FluxA
        FluxB = problem.FluxB
        if self.corrector is not None:
            assert problem.FluxCorr is not None
            FluxC = problem.FluxCorr

        h = (tF - t0) / float(N-1)
        eps = problem.epsilon
        state = x0.copy()
        yield state.copy()

        if self.corrector is not None:
            corr_c = self.corrector.c * (eps**self.corrector.eps_ord) * (h**self.corrector.h_ord)

        #yield state
        for i in range(1, N):
            if self.corrector is not None:
                FluxC(-h * corr_c, state, problem.params)

            # c1 d1 ..... c_{steps} d_{steps}
            for k in range(len(self.cs)):
                FluxA(h*self.cs[k], state, problem.params)
                FluxB(eps*h*self.ds[k], state, problem.params)

            # Follow A flux for c_{step}
            FluxA(h*self.cs[-1], state, problem.params)

            # d{steps-1} c{steps-1} ..... d1 c1
            for k in range(len(self.ds)-2, -1, -1):
                FluxB(eps*h*self.ds[k], state, problem.params)
                FluxA(h*self.cs[k], state, problem.params)

            if self.corrector is not None:
                FluxC(-h * corr_c, state, problem.params)

            yield state.copy()

    def Integrate(self, problem, t0, tF, x0, N):
        FluxA = problem.FluxA
        FluxB = problem.FluxB
        if self.corrector is not None:
            assert problem.FluxCorr is not None
            FluxC = problem.FluxCorr

        h = (tF - t0) / float(N-1)
        eps = problem.epsilon

        hist = np.zeros((N, len(x0)), dtype=x0.dtype)
        hist[0,:] = x0

        if self.corrector is not None:
            corr_c = self.corrector.c * (eps**self.corrector.eps_ord) * (h**self.corrector.h_ord)

        for i in range(1, N):
            state = hist[i-1,:].copy()
            if self.corrector is not None:
                FluxC(-h * corr_c, state, problem.params)

            # c1 d1 ..... c_{steps} d_{steps}
            for k in range(len(self.cs)):
                FluxA(h*self.cs[k], state, problem.params)
                FluxB(eps*h*self.ds[k], state, problem.params)

            # Follow A flux for c_{step}
            FluxA(h*self.cs[-1], state, problem.params)

            # d{steps-1} c{steps-1} ..... d1 c1
            for k in range(len(self.ds)-2, -1, -1):
                FluxB(eps*h*self.ds[k], state, problem.params)
                FluxA(h*self.cs[k], state, problem.params)

            if self.corrector is not None:
                FluxC(-h * corr_c, state, problem.params)

            hist[i,:] = state
            #print(i)

        return hist



cs_leap = np.array([0.5])
ds_leap = np.array([1.0])
LeapFrog = SplittingMethod(cs_leap, ds_leap)

cs_saba3 = np.array([(5. - np.sqrt(15.)) / 10. , np.sqrt(15.) / 10.])
ds_saba3 = np.array([5./18., 4./9.])
SABA3 = SplittingMethod(cs_saba3, ds_saba3)

corr_saba3 = Corrector((54. - 13. * np.sqrt(15.)) / (2. * 648.), 2, 2)
SABA3_corr = SplittingMethod(cs_saba3, ds_saba3, corr_saba3)
