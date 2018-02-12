#include <iostream>
#include <string>
#include <fstream>
#include <queue>
#include <algorithm>
#include "dem.h"
#include "Node.h"
#include "utils.h"
#include <time.h>
#include <list>
#include <stack>
#include <vector>
using namespace std;

typedef std::vector<Node> NodeVector;
typedef std::priority_queue<Node, NodeVector, Node::Greater> PriorityQueue;

void InitPriorityQue(CDEM& dem, Flag& flag,  PriorityQueue& priorityQueue)
{
	int width=dem.Get_NX();
	int height=dem.Get_NY();
	Node tmpNode;
	int iRow, iCol,row,col;

	queue<Node> depressionQue;

	// push border cells into the PQ
	for (row = 0; row < height; row++)
	{
		for (col = 0; col < width; col++)
		{
			if (flag.IsProcessedDirect(row,col)) continue;

			if (dem.is_NoData(row, col)) {
				flag.SetFlag(row,col);
				for (int i = 0; i < 8; i++)
				{
					iRow = Get_rowTo(i, row);
					iCol = Get_colTo(i, col);
					if (flag.IsProcessed(iRow,iCol)) continue;
					if (!dem.is_NoData(iRow, iCol))
					{
						tmpNode.row = iRow;
						tmpNode.col = iCol;
						tmpNode.spill = dem.asFloat(iRow, iCol);
						priorityQueue.push(tmpNode);
						flag.SetFlag(iRow,iCol);
					}
				}
			}
			else
			{
				if (row==0 || row==height-1 || col==0 || col==width-1){
					//on the DEM border
					tmpNode.row = row;
					tmpNode.col = col;
					tmpNode.spill = dem.asFloat(row, col);
					priorityQueue.push(tmpNode);
					flag.SetFlag(row,col);					
				}
			}
		}
	}
}
void ProcessTraceQue(CDEM& dem,Flag& flag,queue<Node>& traceQueue, PriorityQueue& priorityQueue) 
{
	bool HaveSpillPathOrLowerSpillOutlet;
	int i,iRow,iCol;
	int k,kRow,kCol;
	int noderow,nodecol;
	Node N,node;
	queue<Node> potentialQueue;
	int indexThreshold=2;  //index threshold, default to 2
	while (!traceQueue.empty())
	{
		node = traceQueue.front();
		traceQueue.pop();
		noderow=node.row;
		nodecol=node.col;
		bool Mask[5][5]={{false},{false},{false},{false},{false}};
		for (i = 0; i < 8; i++){
			iRow = Get_rowTo(i,noderow);
			iCol = Get_colTo(i,nodecol);
			if(flag.IsProcessedDirect(iRow,iCol)) continue;
			if (dem.asFloat(iRow,iCol)>node.spill){
				N.col = iCol;
				N.row = iRow;
				N.spill = dem.asFloat(iRow,iCol);
				traceQueue.push(N);
				flag.SetFlag(iRow,iCol);
			}
			else{
				//initialize all masks as false		
				HaveSpillPathOrLowerSpillOutlet=false; //whether cell i has a spill path or a lower spill outlet than node if i is a depression cell
				for(k = 0; k < 8; k++){
					kRow = Get_rowTo(k,iRow);
					kCol = Get_colTo(k,iCol);
					if((Mask[kRow-noderow+2][kCol-nodecol+2]) ||
						(flag.IsProcessedDirect(kRow,kCol)&&dem.asFloat(kRow,kCol)<node.spill)
						)
					{
						Mask[iRow-noderow+2][iCol-nodecol+2]=true;
						HaveSpillPathOrLowerSpillOutlet=true;
						break;
					}
				}
				if(!HaveSpillPathOrLowerSpillOutlet){
					if (i<indexThreshold) potentialQueue.push(node);
					else
						priorityQueue.push(node);
					break; // make sure node is not pushed twice into PQ
				}
			}
		}//end of for loop
	}

	while (!potentialQueue.empty())
	{
		node = potentialQueue.front();
		potentialQueue.pop();
		noderow=node.row;
		nodecol=node.col;

		//first case
		for (i = 0; i < 8; i++)
		{
			iRow = Get_rowTo(i,noderow);
			iCol = Get_colTo(i,nodecol);
			if(flag.IsProcessedDirect(iRow,iCol)) continue;
			else {
				priorityQueue.push(node);
				break;
			}
		}		
	}
}

void ProcessPit(CDEM& dem, Flag& flag, queue<Node>& depressionQue,
	queue<Node>& traceQueue,PriorityQueue& priorityQueue)
{
	int iRow, iCol,i;
	float iSpill;
	Node N;
	Node node;
	int width=dem.Get_NX();
	int height=dem.Get_NY();
	while (!depressionQue.empty())
	{
		node= depressionQue.front();
		depressionQue.pop();
		for (i = 0; i < 8; i++)
		{
			iRow = Get_rowTo(i, node.row);
			iCol = Get_colTo(i,  node.col);
			if (flag.IsProcessedDirect(iRow,iCol)) continue;		
			iSpill = dem.asFloat(iRow, iCol);
			if (iSpill > node.spill)     
			{ //slope cell
				N.row = iRow;
				N.col = iCol;
				N.spill = iSpill;				
				flag.SetFlag(iRow,iCol);
				traceQueue.push(N);
				continue;
			}
			//depression cell
			flag.SetFlag(iRow,iCol);
			dem.Set_Value(iRow, iCol, node.spill);
			N.row = iRow;
			N.col = iCol;
			N.spill = node.spill;
			depressionQue.push(N);
		}
	}
}
void fillDEM(char* inputFile, char* outputFilledPath)
{
	queue<Node> traceQueue;
	queue<Node> depressionQue;
	//read float-type DEM
	CDEM dem;
	double geoTransformArgs[6];
	std::cout<<"Reading input tiff file..."<<endl;
	if (!readTIFF(inputFile, GDALDataType::GDT_Float32, dem, geoTransformArgs)){
		printf("Error occurred while reading GeoTIFF file!\n");
		return;
	}	
	std::cout<<"Finish reading file"<<endl;
	time_t timeStart, timeEnd;
	int width = dem.Get_NX();
	int height = dem.Get_NY();
	std::cout<<"Using our proposed variant to fill DEM"<<endl;
	timeStart = time(NULL);
	Flag flag;
	if (!flag.Init(width,height)) {
		printf("Failed to allocate memory!\n");
		return;
	}
	PriorityQueue priorityQueue;
	int iRow, iCol, row,col;
	float iSpill,spill;

	int numberofall=0;
	int numberofright=0;

	InitPriorityQue(dem,flag,priorityQueue); 
	while (!priorityQueue.empty())
	{
		Node tmpNode = priorityQueue.top();
		priorityQueue.pop();
		row = tmpNode.row;
		col = tmpNode.col;
		spill = tmpNode.spill;

		for (int i = 0; i < 8; i++)
		{
			iRow = Get_rowTo(i, row);
			iCol = Get_colTo(i, col);

			if (flag.IsProcessed(iRow,iCol)) continue;
			iSpill = dem.asFloat(iRow, iCol);
			if (iSpill <= spill){
				//depression cell
				dem.Set_Value(iRow, iCol, spill);
				flag.SetFlag(iRow,iCol);
				tmpNode.row = iRow;
				tmpNode.col = iCol;
				tmpNode.spill = spill;
				depressionQue.push(tmpNode);
				ProcessPit(dem,flag,depressionQue,traceQueue,priorityQueue);
			}
			else
			{
				//slope cell
				flag.SetFlag(iRow,iCol);
				tmpNode.row = iRow;
				tmpNode.col = iCol;
				tmpNode.spill = iSpill;
				traceQueue.push(tmpNode);
			}			
			ProcessTraceQue(dem,flag,traceQueue,priorityQueue); 
		}
	}
	timeEnd = time(NULL);
	double consumeTime = difftime(timeEnd, timeStart);
	std::cout<<"Time used:"<<consumeTime<<" seconds"<<endl;

	double min, max, mean, stdDev;
	calculateStatistics(dem, &min, &max, &mean, &stdDev);
	CreateGeoTIFF(outputFilledPath, dem.Get_NY(), dem.Get_NX(), 
		(void *)dem.getDEMdata(),GDALDataType::GDT_Float32, geoTransformArgs,
		&min, &max, &mean, &stdDev, -9999);
	return;
}