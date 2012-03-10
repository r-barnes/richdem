/*
//An adaptation of Paul Heckbert's "A Seed Fill Algorithm" from "Graphics Gems"
#define bucket_push(y,xl,xr,dy) if(y+(dy)>=0 && y+(dy)<height) segments.push(segment(y,xl,xr,dy))
#define bucket_push_hard(y,xl,xr,dy) segments.push(segment(y,xl,xr,dy))

typedef struct segment_type{
	int y, xl, xr, dy;
	segment_type(int y, int xl, int xr, int dy) : y(y+dy), xl(xl), xr(xr), dy(dy) {}
} segment;

void label_this(int x, int y, const int ov, const int nv, int_2d &groups, const char_2d &flowdirs){
	std::stack<segment> segments;
	int l, x1, x2, dy;
	const int width=flowdirs.width();
	const int height=flowdirs.height();

    if (ov==nv || x<0 || x>=width || y<0 || y>=height) return;
	if (groups(x,y)!=ov || groups(x,y)==nv || flowdirs(x,y)!=NO_FLOW) return;
    bucket_push(y,x,x,1);			// needed in some cases
    bucket_push(y+1,x,x,-1);		// seed segment (popped 1st)

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
	
		for (x=x1; x>=0 && groups(x,y)==ov && flowdirs(x,y)==NO_FLOW; x--)
			groups(x, y)=nv;
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
			for (; x<width && groups(x, y)==ov && flowdirs(x,y)==NO_FLOW; x++)
				groups(x, y)=nv;
			bucket_push(y, l, x-1, dy);
			if (x>x2+1)						// leak on right?
				bucket_push(y, x2+1, x-1, -dy);

			//Diagonal Leaks?
			if(x<width){
				if(y+dy>=0 && y+dy<height && groups(x-1,y+dy)!=ov)
					bucket_push_hard(y,x,x,dy);
				if(y-dy>=0 && y-dy<height && groups(x-1,y-dy)!=ov)
					bucket_push_hard(y,x,x,-dy);
			}
			if(x-1==width-1){
				if(y+dy>=0 && y+dy<height && groups(x-1,y+dy)!=ov)
					bucket_push_hard(y,x-2,x-2,dy);
				if(y-dy>=0 && y-dy<height && groups(x-1,y-dy)!=ov)
					bucket_push_hard(y,x-2,x-2,-dy);
			}

skip:		for (x++; x<=x2 && groups(x, y)!=ov && flowdirs(x,y)==NO_FLOW; x++);
			l = x;
		} while (x<=x2);
    }
}*/
