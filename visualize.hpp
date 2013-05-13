#ifndef _visualize_included
#define _visualize_included

#include "CImg.hpp"
#include "data_structures.hpp"
#include "interface.hpp"
#include <iostream>
#include <string>
#include <limits>
using namespace cimg_library;

static const unsigned char white[]={255,255,255}, red[]={255,0,0}, green[]={0,255,0};

template <class T>
void make_image(CImg<unsigned char> &pout, const array2d<T> &data, bool do_highlight, T highlight, T zmin, T zmax, const unsigned char grayscale[]){
  diagnostic("Drawing image...");
  int tc;
  for(int x=0;x<data.width();x++)
    for(int y=0;y<data.height();y++){
//      std::cerr<<x<<" "<<y<<" "<<data(x,y)<<std::endl;
      if(data(x,y)==data.no_data)
        pout.draw_point(x,y,red);
      else if(do_highlight && data(x,y)==highlight)
        pout.draw_point(x,y,green);
      else if (data(x,y)<zmin)
        pout.draw_point(x,y,&grayscale[0]);
      else if (data(x,y)>zmax)
        pout.draw_point(x,y,&grayscale[255*3]);
      else {
        tc=(int)((data(x,y)-zmin)*255*3/(zmax-zmin));
        pout.draw_point(x,y,&grayscale[tc]);
      }
    }
  diagnostic("succeeded.\n");
}

template<class T>
void visualize(const array2d<T> &data, bool do_highlight, T highlight, const char title[]){
  CImg<unsigned char> pout(data.width(),data.height(), 1, 3, 0); //Image only uses top 1024 pixels. Text uses the rest.
  unsigned char grayscale[256*3];

  diagnostic("Constructing grayscale table...");
  #pragma omp parallel for
  for(int i=0;i<256;i++){
    grayscale[i*3+0]=i;
    grayscale[i*3+1]=i;
    grayscale[i*3+2]=i;
  }
  diagnostic("succeeded.\n");

  diagnostic("Searching for maximum and minimum elevations...");
  T zmin=std::numeric_limits<T>::max();
  T zmax=std::numeric_limits<T>::min();
  for(int x=0;x<data.width();x++)
    for(int y=0;y<data.height();y++)
      if(data(x,y)==data.no_data || data(x,y)==highlight)
        continue;
      else if(data(x,y)<zmin)
        zmin=data(x,y);
      else if (data(x,y)>zmax)
        zmax=data(x,y);
  diagnostic("succeeded.\n");

  diagnostic("Displaying image...\n");
  CImgDisplay main_disp(1000,500,title);
  int factor = 400, x = factor, y = factor,x0=x-factor,y0=y-factor,x1=x+factor,y1=y+factor;
  float zmin_scale=1,zmax_scale=1;
  int mx=-1,my=0;
  bool redraw = true, zscale = true;
  while (!main_disp.is_closed()) {
    if (zscale){
      diagnostic_arg("Zmax is now at %.2f%% (%f)\n",zmax_scale*100.0,zmax*zmax_scale);
      make_image(pout, data, do_highlight, highlight, (T) (zmin*zmin_scale), (T)(zmax*zmax_scale), grayscale);
      zscale=false;
    }
    if (redraw) {
      x0=x-factor;
      y0=y-factor;
      x1=x+factor;
      y1=y+factor;

      CImg<unsigned char> visu = pout.get_crop(x0,y0,x1,y1);
      pout.get_crop(x0,y0,x1,y1).resize(main_disp).draw_text(2,2,"Bob",white,0,1,13,x,y).display(main_disp);
      redraw=false;
    }
    if (main_disp.is_keyI()) { factor = (int)(factor/1.2f); if (factor<3) factor = 3; redraw = true; }
    else if (main_disp.is_keyO()) { factor = (int)(factor*1.2f); if (factor>800) factor = 800; redraw = true; }
    else if (main_disp.is_keyARROWUP()) { y-=factor*.1; redraw=true; }
    else if (main_disp.is_keyARROWDOWN()) { y+=factor*.1; redraw=true; }
    else if (main_disp.is_keyARROWLEFT()) { x-=factor*.1; redraw=true; }
    else if (main_disp.is_keyARROWRIGHT()) { x+=factor*.1; redraw=true; }
    else if (main_disp.is_keyZ()) {
      zmin_scale-=0.1;
      zscale=true;
      if(zmin_scale<1)
        zmin_scale=1;
    } else if (main_disp.is_keyX()) {
      if(! (zmin*(zmin_scale+0.1)>zmax) ){
        zmin_scale+=0.1;
        zscale=true;
      }
    } else if (main_disp.is_keyM()) {
      zmax_scale+=0.1;
      zscale=true;
      if(zmax_scale>1)
        zmax_scale=1;
    } else if (main_disp.is_keyN()) {
      if(! (zmax*(zmax_scale-0.05)<zmin) ){
        zmax_scale-=0.05;
        zscale=true;
      }
    } else if (main_disp.is_keyB()){
      std::cout<<"Zmax_scale (0-100%): ";
      std::cin>>zmax_scale;
      if(zmax_scale>1)
        zmax_scale=1;
      if(zmax*zmax_scale<zmin)
        zmax_scale=1;
      zscale=true;
    } else if (main_disp.is_keyS()) {
      std::string fname;
      std::cout<<"Filename for saved image: ";
      std::cin>>fname;
      pout.save_png(fname.c_str());
      std::cout<<"Saved image."<<std::endl;
    } else if (main_disp.is_resized()) { main_disp.resize(); redraw = true; }
    else if (main_disp.button() && main_disp.mouse_x()>=0) {
      mx = main_disp.mouse_x()/((double)main_disp.width())*(x1-x0)+x0;
      my = main_disp.mouse_y()/((double)main_disp.height())*(y1-y0)+y0;
      diagnostic_arg("Clicked (%d,%d).\n",mx,my);
    } else if (main_disp.is_keyE()) //Exit program
      break;
    main_disp.wait();
  }
  diagnostic("succeeded.\n");

  return;
}

#endif
