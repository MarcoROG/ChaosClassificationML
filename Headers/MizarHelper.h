#pragma once

#include <stdio.h>
extern "C"{
    #include "mizar.h"
}
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

//Graphical area
#define GXMAX 1023
#define GYMAX 1023
#define SHIFT 100

using namespace std;

/** 
 * Provides helper functions for the Mizar graphic library
 *
 */
class MizarHelper
{
public:
    /**
     * Initializes the graphical environment.
     *
     * @param min_x     Minimum value of x axis
     * @param min_y     Minimum value of y axis
     * @param max_x     Maximum value of x axis
     * @param max_y     Maximum value of y axis
     * @param px_w      Pixel width of the window.        Default 1023
     * @param px_h      Pixel height of the window.       Default 1023
     * @param marg      Pixel margin around the drawing.  Default 100
     */
    MizarHelper(double min_x, double min_y, double max_x, double max_y, int px_w = GXMAX, int px_h = GYMAX, int marg = SHIFT)
        : px_width(px_w), px_height(px_h), ua_min_x(min_x), ua_min_y(min_y), ua_max_x(max_x), ua_max_y(max_y), margin(marg)
    {
        // Initialize Mizar
        inizio(0); // 0 = how the output works
        
        defpag(2*margin + px_width, 2*margin + px_width); // Page size
        grarea(margin, margin + px_width, margin, margin + px_width); // Graphical area
        utarea(ua_min_x, ua_max_x, ua_min_y, ua_max_y); // User area
    }

    /**
     * Resizes the user coordinates range
     *
     * @param min_x     Minimum value of x axis
     * @param min_y     Minimum value of y axis
     * @param max_x     Maximum value of x axis
     * @param max_y     Maximum value of y axis
     */
    void Resize(double min_x, double max_x, double min_y, double max_y)
    {
        ua_min_x = min_x;
        ua_max_x = max_x;
        ua_min_y = min_y;
        ua_max_y = max_y;

        utarea(ua_min_x, ua_max_x, ua_min_y, ua_max_y); // User area
    }

    /**
     * Creaes a new plotting page with the same dimensions as the previous
     */
    void NewPlot()
    {
        pagina();
    }

    /**
     * Creates a new plotting page with the given dimensions
     *
     * @param min_x     Minimum value of x axis
     * @param min_y     Minimum value of y axis
     * @param max_x     Maximum value of x axis
     * @param max_y     Maximum value of y axis
     */
    void NewPlot(double min_x, double max_x, double min_y, double max_y)
    {
        this->NewPlot();        
        this->Resize(min_x, max_x, min_y, max_y);
    }

    /**
     * Adds linear axes to the current graph
     * 
     * @param lab_x Label for x axis
     * @param lab_y Label for y axis
     */
    void LinearAxes(string lab_x, string lab_y)
    {
        coord(1);
        quadro((char*)lab_x.c_str(), lab_x.length(), (char*) lab_y.c_str(), lab_y.length());
    }

    /**
     * Adds logarithmic axes to the current graph
     * 
     * @param lab_x Label for x axis
     * @param lab_y Label for y axis
     * @param size  Text dimension. Default is 1.0
     */
    void LogLogAxes(string lab_x, string lab_y, double size = 1.0)
    {
        int latoX[2] = {-1, 1};
        int latoY[2] = {1, -1};

        coord(4);
        asslog(ua_min_x, ua_min_y, ua_max_x, ua_min_y, ua_min_x, ua_max_x, latoX, size, (char*)lab_x.c_str(), lab_x.length());
        asslog(ua_min_x, ua_min_y, ua_min_x, ua_max_y, ua_min_y, ua_max_y, latoY, size, (char*)lab_y.c_str(), lab_y.length());
    }
    
    /**
     * Adds a symbol to the specified coordinates
     *
     * @see http://www.mat.unimi.it/users/sansotte/metmod/lezioni/grafico-crop.pdf for symbols table
     * @param symbol    Octal ID of the desired symbol
     * @param x         x position of symbol
     * @param y         y position of symbol
     * @param size      Symbol dimension. Default is 1.0
     * @param angle     Symbol rotation in radians. Default is 0.0
     */
    void PutSymbol(int symbol, double x, double y, double size = 1.0, double angle = 0.0)
    {
        defsim(symbol, size, angle);
        utsim(x, y);
        this->CursorTo(ua_max_x, ua_max_y);
    }

    /**
     * Adds a string of text to the plot
     *
     * @param text  The text to be displayed
     * @param x     x position of text (start)
     * @param y     y position of text (bottom)
     * @param size  Text dimension. Default is 1.0
     * @param angle Text rotation in radians. Default is 0.0
     */
    void PutText(string text, double x, double y, double size = 1.0, double angle = 0.0)
    {
        utmov(x, y);  
        stralf((char*)text.c_str(), text.length(), size, angle);  
    }

    /**
     * Plots the two sequences of data, assuming they have the same length
     *
     * @param xs    X coordinates to be plotted
     * @param ys    Y coordinates to be plotted
     */
    void Plot(vector<double> xs, vector<double> ys)
    {
        this->CursorTo(xs[0], ys[0]);
        
        for (int i = 1; i < xs.size(); i++)
        {
            this->LineTo(xs[i], ys[i]);
        }
    }

    /**
     * Plots two selected sequences from the simulation data. First axis is time, the second is dimension
     * 
     * @param history   Simulation data
     * @param id_x      Dimension index to be plotted on X axis
     * @param id_y      Dimension index to be plotted on Y axis
     */
    void PlotHistory(vector<vector<double>>& history, int id_x, int id_y)
    {
        this->CursorTo(history[0][id_x], history[0][id_y]);

        for (int i = 1; i < history.size(); i++)
        {
            this->LineTo(history[i][id_x], history[i][id_y]);
        }
    }

    /**
     * Moves the cursor to the desired location
     * 
     * @param x     X coordinate of desired location 
     * @param y     Y coordinate of desired location
     */
    void CursorTo(double x, double y)
    {
        utmov(x,y);
    }

    /**
     * Draws a line between cursor and desired location
     * 
     * @param x     X coordinate of desired location 
     * @param y     Y coordinate of desired location
     */
    void LineTo(double x, double y)
    {
        utseg(x,y);
    }

    /**
     *  Reads either coordinates or a character from user input
     *
     * @param   x   The x coordinate of the user click
     * @param   y   The y coordinate of the user click
     * @param   c   The character the user pressed 
     */
    void UserInput(double* x, double* y, char* c)
    {
        utcur(x,y,c);
    }

    /**
     * Waits for user input before stopping
     */
    void Stop()
    {
        fine();
    }

    ~MizarHelper()
    {
        fine();
    }

    /**
     * Truncates the string representation of a double to a given number of decimals
     *
     * @param number    The number to be represented
     * @param decimals  Desired number of decimals
     * @return          String representation of said number with given number of decimals
     */
    string GetTrunc(double number, int decimals)
    {
        ostringstream out;
        out.precision(decimals);
        out << std::fixed << number;
        return out.str();
    }

private:
    int px_width, px_height, margin;
    double ua_min_x, ua_max_x;
    double ua_min_y, ua_max_y;
};
