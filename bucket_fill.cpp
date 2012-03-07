#include <stack>
#include "data_structures.h"
#include "utility.h"
#include <cstdlib>

typedef struct segment_type{
	int y, xl, xr, dy;
	segment_type(int y, int xl, int xr, int dy) : y(y+dy), xl(xl), xr(xr), dy(dy) {}
} segment;

void bucket_fill(int x, int y, int oldv, int nv, char_2d &data){
	std::stack<segment> segments;
	int l, x1, x2, dy;
    int ov;	/* old pixel value */

    ov = data(x, y);		/* read pv at seed point */
    if (ov==nv || x<0 || x>=data.width() || y<0 || y>=data.height()) return;
    segments.push(segment(y, x, x, 1));			/* needed in some cases */
    segments.push(segment(y+1, x, x, -1));		/* seed segment (popped 1st) */

    while (segments.size()>0) {
	/* pop segment off stack and fill a neighboring scan line */
	segment c=segments.top();
	segments.pop();
	y=c.y;
	x1=c.xl;
	x2=c.xr;
	dy=c.dy;
	/*
	 * segment of scan line y-dy for x1<=x<=x2 was previously filled,
	 * now explore adjacent pixels in scan line y
	 */
	for (x=x1; x>=0 && data(x, y)==ov; x--)
	    data(x, y)=nv;
	if (x>=x1) goto skip;
	l = x+1;
	if (l<x1) segments.push(segment(y, l, x1-1, -dy));		/* leak on left? */
	x = x1+1;
	do {
	    for (; x<data.width() && data(x, y)==ov; x++)
		data(x, y)=nv;
	    segments.push(segment(y, l, x-1, dy));
	    if (x>x2+1) segments.push(segment(y, x2+1, x-1, -dy));	/* leak on right? */
skip:	    for (x++; x<=x2 && data(x, y)!=ov; x++);
	    l = x;
	} while (x<=x2);
    }
}

int printdata(const char_2d &data){
	for(int y=0;y<data.height();y++){
		for(int x=0;x<data.width();x++)
			printf("%d ",data(x,y));
		printf("\n");
	}
}

int random_array(int width, int height, char_2d &data){
	data.resize(width,height,false);
	for(int x=0;x<width;x++)
	for(int y=0;y<height;y++)
		data(x,y)=rand()%2;
}

int main(){
	char_2d data;
	random_array(30,30,data);
	data(15,15)=0;
	printdata(data);
	bucket_fill(15,15,0,3,data);
	printdata(data);
}
