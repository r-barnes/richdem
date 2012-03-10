#include <stdio.h>

int main(){
	int i;
	float f;
	for(i=0,f=0;;i++,f++)
		if(i!=f)
			printf("%d %f\n",i,f);
}
