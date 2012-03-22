#ifndef _visualize_included
#define _visualize_included

#include "CImg.h"
#include "data_structures.h"
using namespace cimg_library;

template<class T, class U>
void visualize(const array2d<T> &data, bool do_highlight, U highlight){
	CImg<unsigned char> pout(data.width(),data.height(), 1, 3, 0); //Image only uses top 1024 pixels. Text uses the rest.
	const unsigned char white[]={255,255,255}, red[]={255,0,0}, green[]={0,255,0};
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
	float zmin=9e99,zmax=-9e99;
	for(int x=0;x<data.width();x++)
		for(int y=0;y<data.height();y++)
			if(data(x,y)==data.no_data || data(x,y)==highlight)
				continue;
			else if(data(x,y)<zmin)
				zmin=data(x,y);
			else if (data(x,y)>zmax)
				zmax=data(x,y);
	diagnostic("succeeded.\n");

	diagnostic("Drawing image...");
	int tc;
	for(int x=0;x<data.width();x++)
		for(int y=0;y<data.height();y++){
			if(data(x,y)==data.no_data)
				pout.draw_point(x,y,green);
			else if(do_highlight && data(x,y)==highlight)
				pout.draw_point(x,y,red);
			else{
				tc=(int)((data(x,y)-zmin)*255*3/(zmax-zmin));
				pout.draw_point(x,y,&grayscale[tc]);
			}
		}
	diagnostic("succeeded.\n");

	diagnostic("Displaying image...");
	CImgDisplay main_disp(1000,500,"RichDEM Visualization");
	int factor = 400, x = factor, y = factor, mx=-1,my=0,x0=x-factor,y0=y-factor,x1=x+factor,y1=y+factor;
	bool redraw = true;
	while (!main_disp.is_closed()) {
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
		if (main_disp.is_keyO()) { factor = (int)(factor*1.2f); if (factor>800) factor = 800; redraw = true; }
		if (main_disp.is_keyARROWUP()) { y-=factor*.1; redraw=true; }
		if (main_disp.is_keyARROWDOWN()) { y+=factor*.1; redraw=true; }
		if (main_disp.is_keyARROWLEFT()) { x-=factor*.1; redraw=true; }
		if (main_disp.is_keyARROWRIGHT()) { x+=factor*.1; redraw=true; }
		if (main_disp.is_resized()) { main_disp.resize(); redraw = true; }
		if (main_disp.button() && main_disp.mouse_x()>=0) {
			mx = main_disp.mouse_x()/((double)main_disp.width())*(x1-x0)+x0;
			my = main_disp.mouse_y()/((double)main_disp.height())*(y1-y0)+y0;
		} else if (main_disp.is_keyE()) //Exit program
			break;
		main_disp.wait();
	}
	diagnostic("succeeded.\n");

	return;
}

#endif
