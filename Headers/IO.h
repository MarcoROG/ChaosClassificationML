#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


/**
 * Prints a given matrix to file in text form
 *
 * @param fileName The file path
 */
void SaveMatrix(string fileName, vector<vector<double>>& matrix, int skip)
{
    ofstream file(fileName);
    file << std::fixed << std::setprecision(15);

    int nrows = matrix.size();
    int ncols = matrix[0].size();

    for (int i = 0; i < nrows; i+=skip)
    {
        for (int j = 0; j < ncols; j++)
            file << matrix[i][j] << "\t";

        file << std::endl;
    }

    file.close();
}

/**
 * Prints a given matrix to file in binary form
 *
 * @param fileName The file path
 */
void SaveMatrixBin(string fileName, vector<vector<double>>& matrix, int skip)
{
    ofstream file(fileName, ios::binary | ios::out);
    file << std::fixed << std::setprecision(15);

    int nrows = matrix.size();
    int ncols = matrix[0].size();

    int nlines = nrows / skip;

    file.write(reinterpret_cast<const char *>(&nlines), sizeof(nlines));
    file.write(reinterpret_cast<const char *>(&ncols), sizeof(ncols));

    for (int i = 0; i < nrows; i+=skip)
    {
        // Store its contents
        file.write(reinterpret_cast<const char *>(&matrix[i][0]), matrix[i].size()*sizeof(double));
    }

    file.close();
}

inline void SaveRowBin(ofstream& file, double* row, int n)
{  
    file.write(reinterpret_cast<const char *>(row), n*sizeof(double));
}
