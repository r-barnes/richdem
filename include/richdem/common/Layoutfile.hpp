/**
  @file
  @brief Defines classes used for reading and writing tiled datasets.

  A layout file is a text file with the format:

      file1.tif, file2.tif, file3.tif,
      file4.tif, file5.tif, file6.tif, file7.tif
               , file8.tif,          ,

  where each of fileX.tif is a tile of the larger DEM collectively described by
  all of the files. All of fileX.tif must have the same shape; the layout file
  specifies how fileX.tif are arranged in relation to each other in space.
  Blanks between commas indicate that there is no tile there: the algorithm will
  treat such gaps as places to route flow towards (as if they are oceans). Note
  that the files need not have TIF format: they can be of any type which GDAL
  can read. Paths to fileX.tif are taken to be relative to the layout file.

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_layoutfile_hpp_
#define _richdem_layoutfile_hpp_

#include "richdem/common/logger.hpp"
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <stdexcept>

namespace richdem {

//Define operating system appropriate directory separators
#if defined(__unix__) || defined(__linux__) || defined(__APPLE__)
  #define RICHDEM_SLASH_CHAR "/"
#elif defined(__WIN32__)
  #define RICHDEM_SLASH_CHAR "\\"
#endif

///Eliminate spaces from the beginning and end of str
static std::string trimStr(std::string const& str){
  if(str.empty())
      return str;

  std::size_t firstScan = str.find_first_not_of(' ');
  std::size_t first     = firstScan == std::string::npos ? str.length() : firstScan;
  std::size_t last      = str.find_last_not_of(' ');
  return str.substr(first, last-first+1);
}

///Get only the filename without its extension. That is, convert
///"path/to/file.ext" to "file"
static std::string GetBaseName(std::string filename){
  auto last_slash  = filename.find_last_of(RICHDEM_SLASH_CHAR);
  auto last_period = filename.find_last_of(".");
  if(last_period!=std::string::npos)
    filename.replace(last_period, std::string::npos, "");
  if(last_slash!=std::string::npos)
    filename.replace(0,last_slash+1,"");
  return filename;
}

///@brief Used for reading a layoutfile describing a tiled dataset.
///
///The class acts as a generator. The layoutfile is read on construction and its
///contents retrieved with next(). The Layoutfile specification can be found in
///Layoutfile.hpp.
class LayoutfileReader {
 private:
  ///Stores the grid of filenames extracted from the layoutfile
  std::vector< std::vector< std::string > > fgrid;

  int gridy        = -1;    //Gets incremented to 0 right away
  int gridx        = -2;    //Gets incremented to -1 right away, which then indicates a new row, the first row
  int new_row      = false; //If true, then the current entry is the beginning of a new row. Gets set to true right away.

  std::string filename;     //Of "path/to/file.ext" this is "file.ext"
  std::string basename;     //Of "path/to/file.ext" this is "file"
  std::string path;         //Of "path/to/file.ext" this is "path/to/"
 public:

  ///@brief  Construct a new LayoutfileReader object reading from a given file.
  ///@author Richard Barnes
  ///
  ///@param  layout_filename Layoutfile to read from.
  LayoutfileReader(std::string layout_filename){
    assert(layout_filename.size()>0);
    path = "";
    std::size_t last_slash = layout_filename.find_last_of(RICHDEM_SLASH_CHAR);
    if(last_slash!=std::string::npos)
      path = layout_filename.substr(0,last_slash+1);
    RDLOG_CONFIG<<"Base path for layout-identified files = "<<path;

    //Open file
    std::ifstream fin_layout(layout_filename);

    //Did the file open?
    if(!fin_layout.good())
      throw std::runtime_error("Problem opening layout file '"+layout_filename+"'");

    //Read the entire file
    while(fin_layout){
      //Get a new line from the file.
      std::string row;
      if(!getline(fin_layout, row))
        break;

      //Add another row to the grid to store the data from this line
      fgrid.emplace_back();
      //Put the data somewhere we can process it
      std::istringstream ss(row);

      //While there is unprocessed data
      while(ss){
        //Create temporary item to store incoming datum
        std::string item;
        //Get the next comma-delimited item from the line. If there's data
        //followed by no comma, this returns the data. If there's no data and no
        //further comma, then we'll end with the line one item shorter than it
        //should be. We'll fix this later.
        if (!getline( ss, item, ',' ))
          break;
        //Add the item to the row, trimming off the whitespace
        fgrid.back().push_back(trimStr(item));
      }
    }
    //If we break out of the above loop but haven't reached eof(), then
    //something went wrong.
    if(!fin_layout.eof())
      throw std::runtime_error("Failed to read the entire layout file!");

    //Let's find the longest row
    auto max_row_length = fgrid.front().size();
    for(const auto &row: fgrid)
      max_row_length = std::max(max_row_length,row.size());

    for(auto &row: fgrid)
      if(row.size()==max_row_length-1)   //If the line was one short of max, assume it ends blank
        row.emplace_back();
      else if(row.size()<max_row_length) //The line was more than one short of max: uh oh
        throw std::runtime_error("Not all of the rows in the layout file had the same number of columns!");
  }

  ///@brief Advance the reader to the next layoutfile entry.
  ///@return True if reader advanced successfully; false if not.
  bool next(){
    new_row = false;
    gridx++;

    if(gridx==-1 || gridx==(int)fgrid.front().size()){
      gridy++;
      gridx   = 0;
      new_row = true;
      if(gridy==(int)fgrid.size())
        return false;
    }

    filename = fgrid[gridy][gridx];
    basename = GetBaseName(filename);

    return true;
  }

  ///@return True if the current entry is the beginning of a new row.
  bool newRow() const {
    return new_row;
  }

  ///@return The current entry's filename without the path (e.g. "file.ext").
  const std::string& getFilename() const {
    return filename;
  }

  ///@return The current entry's filename without the path or extension (e.g. "file").
  const std::string& getBasename() const {
    return basename;
  }

  ///@return The current entry's path + filename (e.g. "path/to/file.ext").
  const std::string  getFullPath() const {
    return path + filename;
  }

  ///@brief Return a string representation of the current entry's coordinates.
  ///
  ///A layoutfile is a 2D grid of file names. This method returns the current
  ///entry's position in that grid as `<X>_<Y>`
  ///
  ///@return Current entry's position as a string of the form <X>_<Y>
  const std::string  getGridLocName() const {
    return std::to_string(gridx)+"_"+std::to_string(gridy);
  }

  ///@return Path of layoutfile: of "path/to/layoutfile.layout" returns "path/to/".
  const std::string& getPath() const {
    return path;
  }

  ///@return True if the current entry was a blank
  bool isNullTile() const {
    return filename.size()==0;
  }

  ///@return X-coordinate of the current entry.
  int getX() const {
    return gridx;
  }

  ///@return Y-coordinate of the current entry.
  int getY() const {
    return gridy;
  }
};

///@brief Used for creating a layoutfile describing a tiled dataset.
///
///The class acts as an inverse generator. The layoutfile is created on
///construction and its contents appended to with addEntry(). The Layoutfile
///specification can be found in Layoutfile.hpp.
class LayoutfileWriter {
 private:
  int gridx;           ///Current column being written to
  int gridy;           ///Current row being written to
  std::string path;    ///Path of layoutfile
  std::ofstream flout; ///File output stream for the layout file
 public:

  ///@brief Constructs a new writer object
  ///@param layout_filename Path+Filename of layoutfile to write
  LayoutfileWriter(std::string layout_filename){
    path  = "";
    gridx = 0;
    gridy = 0;
    
    std::size_t last_slash  = layout_filename.find_last_of(RICHDEM_SLASH_CHAR);
    std::size_t last_period = layout_filename.find_last_of(".");
    if(last_slash!=std::string::npos)
      path = layout_filename.substr(0,last_slash+1);
    if(last_period!=std::string::npos)
      layout_filename.replace(last_period+1, std::string::npos, "layout");

    flout.open(layout_filename);
  }

  ~LayoutfileWriter(){
    flout<<std::endl;
  }

  ///@brief Adds a new row to the layoutfile
  void addRow(){
    if(gridx==0 && gridy==0)
      return;
    flout<<std::endl;
    gridy++;
    gridx = 0;
  }

  ///@brief Add a new entry to the layout file.
  ///@param filename File to add. Use `filename=""` to indicate a null tile.
  void addEntry(std::string filename){
    //Get only the filename, not the path to it
    std::size_t last_slash = filename.find_last_of(RICHDEM_SLASH_CHAR);
    if(last_slash!=std::string::npos)
      filename = filename.substr(last_slash+1,std::string::npos);

    if(gridx>0)
      flout<<",";
    flout<<filename;
    gridx++;
  }
};

}

#endif
