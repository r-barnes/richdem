#ifndef _layoutfile_hpp_
#define _layoutfile_hpp_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <stdexcept>

//Define operating system appropriate directory separators
#if defined(__unix__) || defined(__linux__) || defined(__APPLE__)
  #define SLASH_CHAR "/"
#elif defined(__WIN32__)
  #define SLASH_CHAR "\\"
#endif

static std::string trimStr(std::string const& str){
  if(str.empty())
      return str;

  std::size_t firstScan = str.find_first_not_of(' ');
  std::size_t first     = firstScan == std::string::npos ? str.length() : firstScan;
  std::size_t last      = str.find_last_not_of(' ');
  return str.substr(first, last-first+1);
}

static std::string GetBaseName(std::string filename){
  auto last_slash  = filename.find_last_of(SLASH_CHAR);
  auto last_period = filename.find_last_of(".");
  if(last_period!=std::string::npos)
    filename.replace(last_period, std::string::npos, "");
  if(last_slash!=std::string::npos)
    filename.replace(0,last_slash+1,"");
  return filename;
}

class LayoutfileReader {
 private:
  //Stores the grid of filenames
  std::vector< std::vector< std::string > > fgrid;

  int gridy        = -1; //Gets incremented to 0 right away
  int gridx        = -2; //Gets incremented to -1 right away, which then indicates a new row, the first row
  int new_row      = false;

  std::string filename;
  std::string basename;
  std::string path;
 public:
  LayoutfileReader(std::string layout_filename){
    assert(layout_filename.size()>0);
    path = "";
    std::size_t last_slash = layout_filename.find_last_of(SLASH_CHAR);
    if(last_slash!=std::string::npos)
      path = layout_filename.substr(0,last_slash+1);
    std::cerr<<"Base path for layout-identified files: "<<path<<std::endl;

    //Open file
    std::ifstream fin_layout(layout_filename);

    //Did the file open?
    if(!fin_layout.good())
      throw std::runtime_error("Problem opening layout file!");

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

  bool newRow() const {
    return new_row;
  }

  const std::string& getFilename() const {
    return filename;
  }

  const std::string& getBasename() const {
    return basename;
  }

  const std::string  getFullPath() const {
    return path + filename;
  }

  const std::string  getGridLocName() const {
    return std::to_string(gridx)+"_"+std::to_string(gridy);
  }

  const std::string& getPath() const {
    return path;
  }

  bool isNullTile() const {
    return filename.size()==0;
  }

  int getX() const {
    return gridx;
  }

  int getY() const {
    return gridy;
  }
};

class LayoutfileWriter {
 private:
  int gridx;
  int gridy;
  std::string path;
  std::ofstream flout;
 public:
  LayoutfileWriter(std::string layout_filename){
    path  = "";
    gridx = 0;
    gridy = 0;
    
    std::size_t last_slash  = layout_filename.find_last_of(SLASH_CHAR);
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

  void addRow(){
    if(gridx==0 && gridy==0)
      return;
    flout<<std::endl;
    gridy++;
    gridx = 0;
  }

  //Use filename="" to indicate a null tile
  void addEntry(std::string filename){
    //Get only the filename, not the path to it
    std::size_t last_slash = filename.find_last_of(SLASH_CHAR);
    if(last_slash!=std::string::npos)
      filename = filename.substr(last_slash+1,std::string::npos);

    if(gridx>0)
      flout<<",";
    flout<<filename;
    gridx++;
  }
};

#endif