#include <fstream>
#include <map>
#include <sstream>
#include <string>

using namespace std;

typedef map<std::string, double> Parameters;
 

/**
 * Reads the parameters in KEY=VALUE format from a given file
 *
 * @param fileName The file path
 * @return Dictionary<string,double>
 */
Parameters ReadParameters(string fileName)
{
    Parameters params;
    fstream file(fileName);
    string line;
    
    // Read line by line
    while (getline(file, line))
    {
        istringstream is_line(line);
        string key;

        // Read string up to the "="
        if (getline(is_line, key, '='))
        {
            string value;
            // Ignore if it is a comment
            if (key[0] == '#')
                continue;
 
            // Read the value and convert it to double
            if (getline(is_line, value))
            {
                params[key] = atof(value.c_str());
            }
        }
    }
    
    return params;
}
