//Given a flow direction and an area, this tells you what cells receive the area, what proportions of the area, and how much area they got
#include <cmath>
#include <assert.h>
#include <cstdio>
#include <cstdlib>

#define ROUND(A)	floor((A) + 0.5)

//321
//4 0
//567
void where_do_i_flow(float flowdir, int &nhigh, int &nlow){
	//If it is very close to being directed into just one cell
	//then we direct it into just one cell. If we mistakenly direct
	//it into 2 cells, then we may create unresolvable loops in the
	//flow accumulation algorithm, whereas nothing is lost if we
	//leave out one of the two cells (provided a negligible amount
	//of flow is directed to the one we leave out).
	assert(flowdir>=0 && flowdir<=2*M_PI+1e-6);

	float temp=flowdir/(M_PI/4.);

	if(fabs(temp-(int)temp)<1e-6){
		nlow=-1;
		nhigh=(int)ROUND(temp);
	} else {
		nlow=(int)temp;
		nhigh=nlow+1;
	}

	//8 is not technically a direction, but, since things move in a circle,
	//it overlaps with 0. It should _never_ be greater than 8.
	assert(nhigh>=0 && nhigh<=8);
}

//This reacts correctly if the flow direction wedge number exceeds 7.
void area_proportion(float flowdir, int nhigh, int nlow, float &phigh, float &plow){
	if(nlow==-1){
		phigh=1;
		plow=0;
	} else {
		phigh=(nhigh*(M_PI/4.0)-flowdir)/(M_PI/4.0);
		plow=1-phigh;
	}

	assert(phigh+plow==1);	//TODO: This isn't necessarily so in floating-point... or is it?
}

int main(int argc, char **argv){
	if(argc!=3){
		printf("Syntax: PROG flowdir area\n");
		return 0;
	}
	int nhigh, nlow;
	float phigh,plow;
	where_do_i_flow(atof(argv[1]),nhigh,nlow);
	area_proportion(atof(argv[1]),nhigh,nlow,phigh,plow);
	if(nlow==-1)
		printf("Flows to neighbour %d with proportion %f giving %f.\n",nhigh,phigh,atof(argv[2])*phigh);
	else
		printf("Flows to neighbours %d and %d with proportions %f and %f, giving %f and %f.\n",nhigh,nlow,phigh,plow,atof(argv[2])*phigh,atof(argv[2])*plow);

	return 0;
}
