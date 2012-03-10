#include <unistd.h>
#include <stdio.h>

long getTotalSystemMemory()
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
	long pages_free=sysconf(_SC_AVPHYS_PAGES);
    return pages_free * page_size;
}

int main(){
	printf("%ld\n",getTotalSystemMemory());
}
