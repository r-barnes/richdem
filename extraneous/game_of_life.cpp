#include "cpptk.h"
#include <iostream>

using namespace Tk;
using namespace std;

// parameters of the array
int const arrayWidth = 30;
int const arrayHeight = 30;
int squareSize = 10;

// logical array (on/off)
bool cells[arrayWidth][arrayHeight];

// array of canvas elements id
string squares[arrayWidth][arrayHeight];

void setCell(int i, int j, bool state)
{
     ::cells[i][j] = state;
     if (state)
     {
          ".c" << itemconfigure(squares[i][j]) -Tk::fill("red");
     }
     else
     {
          ".c" << itemconfigure(squares[i][j]) -Tk::fill("white");
     }
}

// clears the whole array
void clear()
{
     for (int i = 0; i != arrayWidth; ++i)
     {
          for (int j = 0; j != arrayHeight; ++j)
          {
               setCell(i, j, false);
          }
     }
}

// changes the state of some cell (in response to mouse click)
void click(int x, int y)
{
     // find the logical coordinates
     int i = (x - 1) / squareSize;
     int j = (y - 1) / squareSize;
     
     // toggle the state
     bool newState = !::cells[i][j];
     setCell(i, j, newState);
}

// computes the next generation
void nextStep()
{
     int neighbours[arrayWidth][arrayHeight];

     // initialize the neighbours counters
     for (int i = 1; i != arrayWidth - 1; ++i)
     {
          for (int j = 1; j != arrayHeight - 1; ++j)
          {
               neighbours[i][j] = 0;
          }
     }

     // count the neighbours of each cell
     for (int i = 1; i != arrayWidth - 1; ++i)
     {
          for (int j = 1; j != arrayHeight - 1; ++j)
          {
               if (::cells[i-1][j-1]) ++neighbours[i][j];
               if (::cells[i  ][j-1]) ++neighbours[i][j];
               if (::cells[i+1][j-1]) ++neighbours[i][j];
               if (::cells[i-1][j  ]) ++neighbours[i][j];
               if (::cells[i+1][j  ]) ++neighbours[i][j];
               if (::cells[i-1][j+1]) ++neighbours[i][j];
               if (::cells[i  ][j+1]) ++neighbours[i][j];
               if (::cells[i+1][j+1]) ++neighbours[i][j];
          }
     }
     
     // update the cells (kill or give birth)
     for (int i = 1; i != arrayWidth - 1; ++i)
     {
          for (int j = 1; j != arrayHeight - 1; ++j)
          {
               if (::cells[i][j] == false && neighbours[i][j] == 3)
               {
                    // new cell is born
                    setCell(i, j, true);
               }
               else if (::cells[i][j] == true &&
                    (neighbours[i][j] == 2 || neighbours[i][j] == 3))
               {
                    // remains alive
               }
               else
               {
                    // dies from overcrowding or loneliness or remains dead
                    setCell(i, j, false);
               }
          }
     }
}

int main(int, char *argv[])
{
     try
     {     
          init(argv[0]);
          
          // create the control buttons
          
          frame(".f") -relief(raised) -borderwidth(1);
          button(".f.clear") -text("Clear") -command(::clear);
          button(".f.next") -text("Next") -command(nextStep);
          pack(".f") -side(bottom) -Tk::fill(x);
          pack(".f.clear", ".f.next") -side(Tk::left) -pady(5) -expand(true);
          
          // create the canvas widget
          
          canvas(".c") -background("white")
               -width(squareSize * arrayWidth)
               -height(squareSize * arrayHeight);
          pack(".c") -side(top);
          
          // create and initialize the array of cells
          
          for (int i = 0; i != arrayWidth; ++i)
          {
               for (int j = 0; j != arrayHeight; ++j)
               {
                    ::cells[i][j] = false;
                    
                    Point p1(i * squareSize, j * squareSize);
                    Point p2((i + 1) * squareSize, (j + 1) * squareSize);
                    
                    string squareId(
                         ".c" << create(rectangle, p1, p2)
                         -outline("black") -Tk::fill("white")
                    );
                    squares[i][j] = squareId;
               }
          }
          
          // bind the mouse click so that it is possible
          // to interact with the cells
          
          bind(".c", "<Button-1>", click, event_x, event_y);
          
          wm(resizable, ".", false, false);
          
          runEventLoop();
     }
     catch (exception const &e)
     {
          cerr << "Error: " << e.what() << '\n';
     }
}