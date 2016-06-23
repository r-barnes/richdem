#ifndef _layoutfile_hpp_
#define _layoutfile_hpp_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>

//Define operating system appropriate directory separators
#if defined(__unix__) || defined(__linux__) || defined(__APPLE__)
  #define SLASH_CHAR "/"
#elif defined(__WIN32__)
  #define SLASH_CHAR "\\"
#endif

std::string trimStr(std::string const& str){
  if(str.empty())
      return str;

  std::size_t firstScan = str.find_first_not_of(' ');
  std::size_t first     = firstScan == std::string::npos ? str.length() : firstScan;
  std::size_t last      = str.find_last_not_of(' ');
  return str.substr(first, last-first+1);
}

class LayoutfileReader {
 private:
  std::ifstream     fin_layout;
  std::stringstream cells;
  int gridy     = -1;
  int gridx     = -1;
  int new_row   = false;
  int row_width = -1;
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

    //Read each line of the layout file
    fin_layout.open(layout_filename);

    if(!fin_layout.good())
      throw std::runtime_error("Problem opening layout file!");
  }

  bool next(){
    new_row=false;
    if(cells.peek() == decltype(cells)::traits_type::eof()){
      gridy++;
      new_row = true;

      if(row_width==-1)
        gridx = -1;
      else if(row_width!=gridx)
        throw std::runtime_error("Layout file's rows were not all of the same width!"); //TODO: Print out expected value
      else
        gridx = -1;
      
      std::string line;
      std::getline(fin_layout, line); //Read a line from the layout

      if(line.find_first_not_of("\t\n ")==std::string::npos)
        return false;

      cells = std::stringstream(line);
    }

    std::getline(cells,filename,',');
    gridx++;
    filename = trimStr(filename);

    basename = filename;
    basename.replace(basename.find_last_of("."), std::string::npos, "");
    std::size_t last_slash = basename.find_last_of(SLASH_CHAR);
    if(last_slash!=std::string::npos)
      basename.replace(0,last_slash+1,"");

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

  const std::string& getPath() const {
    return path;
  }

  int getX() const {
    return gridx;
  }

  int getY() const {
    return gridy;
  }
};

#endif