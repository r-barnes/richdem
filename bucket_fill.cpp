#include <stack>
#include "data_structures.h"
#include "utility.h"
#include <cstdlib>

typedef struct segment_type{
	int y, xl, xr, dy;
	segment_type(int y, int xl, int xr, int dy) : y(y+dy), xl(xl), xr(xr), dy(dy) {}
} segment;

#define bucket_push(y,xl,xr,dy) if(y+(dy)>=0 && y+(dy)<height) segments.push(segment(y,xl,xr,dy))
#define bucket_push_hard(y,xl,xr,dy) segments.push(segment(y,xl,xr,dy))

void bucket_fill(int x, int y, int oldv, int nv, char_2d &data){
	std::stack<segment> segments;
	int l, x1, x2, dy;
    int ov;	//old pixel value
	const int width=data.width();
	const int height=data.height();

    ov = data(x, y);		//read pv at seed point
    if (ov==nv || x<0 || x>=width || y<0 || y>=height) return;
    bucket_push(y, x, x, 1);			// needed in some cases
    bucket_push(y+1, x, x, -1);		// seed segment (popped 1st)

    while (segments.size()>0) {
		//pop segment off stack and fill a neighboring scan line
		segment c=segments.top();
		segments.pop();
		y=c.y;
		x1=c.xl;
		x2=c.xr;
		dy=c.dy;
		//segment of scan line y-dy for x1<=x<=x2 was previously filled,
		//now explore adjacent pixels in scan line y
	
		for (x=x1; x>=0 && data(x, y)==ov; x--)
			data(x, y)=nv;
		if (x>=x1){
			if(x-1>=0)							//Diagonal leak?
				bucket_push(y, x-1, x1-1, dy);
			goto skip;
		}
		l = x+1;
		if (l<x1)								// leak on left?
			bucket_push(y, l, x1-1, -dy);

		//Diagonal Leaks?
		if(l-1>=0){
			bucket_push(y, l-1, l-1, dy);
			bucket_push(y, l-1, l-1, -dy);
		}

		x = x1+1;
		do {
			for (; x<width && data(x, y)==ov; x++)
				data(x, y)=nv;
			bucket_push(y, l, x-1, dy);
			if (x>x2+1)						// leak on right?
				bucket_push(y, x2+1, x-1, -dy);

			//Diagonal Leaks?
			if(x<width){
				if(y+dy>=0 && y+dy<height && data(x-1,y+dy)!=ov)
					bucket_push_hard(y,x,x,dy);
				if(y-dy>=0 && y-dy<height && data(x-1,y-dy)!=ov)
					bucket_push_hard(y,x,x,-dy);
			}
			if(x-1==width-1){
				if(y+dy>=0 && y+dy<height && data(x-1,y+dy)!=ov)
					bucket_push_hard(y,x-2,x-2,dy);
				if(y-dy>=0 && y-dy<height && data(x-1,y-dy)!=ov)
					bucket_push_hard(y,x-2,x-2,-dy);
			}

skip:		for (x++; x<=x2 && data(x, y)!=ov; x++);
			l = x;
		} while (x<=x2);
    }
}

typedef struct barnes_segment_type{
	int y, xl, xr, dy, xo;
	barnes_segment_type(int y, int xl, int xr, int dy, int xo) : y(y+dy), xl(xl), xr(xr), dy(dy), xo(xo) {}
} barnes_segment;

//#define push(y,xl,xr,dy,xo) if(y+(dy)>=0 && y+(dy)<data.height()) segments.push(barnes_segment(y,xl,xr,dy,xo))

#define pushd(y,xl,xr,dy) if(IN_GRID(x,y+(dy),data.width(),data.height())) segements.push(segment(y,xl,xr,dy))
/*
void bucket_fill(int seed_x, int seed_y, int oldv, int newv, char_2d &data){
	std::stack<segment> segments;
	int l,x,r;

    if (data(seed_x,seed_y)!=oldv) return;
	if(!IN_GRID(seed_x,seed_y,data.width(),data.height())) return;
    push(seed_y,seed_x,seed_x,1,seed_x);	// needed in some cases
    push(seed_y+1,seed_x,seed_x,-1,seed_x);	// seed segment (popped 1st)

    while (segments.size()>0) {
		//pop segment off stack and fill a neighboring scan line
		barnes_segment c=segments.top();
		segments.pop();
//		fprintf(stderr,"Popped a line %d beginning at %d. Next segments are in %d.\n",c.y,c.xo,c.y+c.dy);
		//segment of scan line y-dy for x1<=x<=x2 was previously filled,
		//now explore adjacent pixels in scan line y
		register bool is_old,is_old_opp;
		int added_to_stack=0;

		x=c.xo;
		for(;x>=0 && data(x,c.y)==oldv;x--){
			if(c.y+c.dy>=0 && c.y+c.dy<data.height()){
				if(is_old && data(x,c.y+c.dy)!=oldv){
					push(c.y,c.xl,c.xr,c.dy,x+1);
					added_to_stack++;
					is_old=false;
				} else if (!is_old && data(x,c.y+c.dy)==oldv)
					is_old=true;
			}
			if(x==c.xl-1)
				is_old_opp=(data(x,c.y-c.dy)==oldv);
			if(x<c.xl && c.y-c.dy>=0 && c.y-c.dy<data.height()){
				if(is_old_opp && data(x,c.y-c.dy)!=oldv){
					push(c.y,c.xl,c.xr,-c.dy,x+1);
					added_to_stack++;
					is_old_opp=false;
				} else if (!is_old_opp && data(x,c.y-c.dy)==oldv)
					is_old_opp=true;
			}
			data(x,c.y)=newv;
		}
		l=x+1;
		if(is_old && c.y+c.dy>=0 && c.y+c.dy<data.height())
			if(data(l,c.y+c.dy)==oldv)
				push(c.y,c.xl,c.xr,c.dy,l);
		if(x<c.xl && is_old_opp && c.y-c.dy>=0 && c.y-c.dy<data.height())
			if(data(l,c.y-c.dy)==oldv)
				push(c.y,c.xl,c.xr,-c.dy,l);


//		fprintf(stderr,"\tLine bounded by [%d,%d].\n",l,r);
	}
}*/

int printdata(const char_2d &data){
	printf(" \t");
	for(int x=0;x<data.width();x++)
		printf("%d",x%10);
	printf("\n");
	for(int y=0;y<data.height();y++){
		printf("%d\t",y%10);
		for(int x=0;x<data.width();x++){
			if(data(x,y)==0)
				printf("\033[31m");
			else if (data(x,y)==1)
				printf("\033[32m");
			else if (data(x,y)==2)
				printf("\033[33m");
			else if (data(x,y)==3)
				printf("\033[34m");
			else if (data(x,y)==4)
				printf("\033[35m");
			printf("%d",data(x,y));
		}
		printf("\n\033[39m");
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
	for(int i=0;i<100;i++){
		random_array(1000,1000,data);
		data(50,50)=0;
		bucket_fill(50,50,0,3,data);
//		printdata(data);
//		printf("==================\n");
	}
}
