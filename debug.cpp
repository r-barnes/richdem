#include <deque>
#include "data_structures.h"
#include "interface.h"

void print_edges(const float_2d &elevations, const std::deque<grid_cell> &low_edges, const std::deque<grid_cell> &high_edges){
	for(int y=0;y<elevations.height();y++){
		for(int x=0;x<elevations.width();x++){
			for(std::deque<grid_cell>::const_iterator i=low_edges.begin();i!=low_edges.end();i++)
				if(x==i->x && y==i->y){
					diagnostic("\033[36m");
					goto printedges_done;
				}
			for(std::deque<grid_cell>::const_iterator i=high_edges.begin();i!=high_edges.end();i++)
				if(x==i->x && y==i->y){
					diagnostic("\033[31m");
					goto printedges_done;
				}
			printedges_done: diagnostic_arg("%2.0f\033[39m ",elevations(x,y));
		}
		diagnostic("\n");
	}
}





template <class T>
void array2d<T>::print_block(std::ostream& out, int minx, int maxx, int miny, int maxy, int precision, std::streamsize swidth){
	out.setf(std::ios::fixed,std::ios::floatfield);
	out<<std::setprecision(precision);
	for(int y=((miny>0)?miny:0);y<=maxy && y<height();y++){
		for(int x=((minx>0)?minx:0);x<=maxx && x<width();x++)
			if(boost::numeric::ublas::matrix<T>::operator()(x,y)==no_data)
				out<<std::setw(swidth)<<"-"<<" ";
			else if(sizeof(T)==1)	//TODO: An ugly way of detecting chars
				out<<std::setw(swidth)<<(int)boost::numeric::ublas::matrix<T>::operator()(x,y)<<" ";
			else
				out<<std::setw(swidth)<<boost::numeric::ublas::matrix<T>::operator()(x,y)<<" ";
		out<<std::endl;
	}
}
