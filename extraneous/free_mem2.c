#include <stdio.h> 
#include <sys/sysinfo.h> 

int main(void) 
{ 
  struct sysinfo myinfo; 
  unsigned long total_bytes; 

  sysinfo(&myinfo); 

  total_bytes = myinfo.mem_unit * myinfo.freeram; 


  printf("total usable main memory is %lu B, %lu MB\n", 
         total_bytes, total_bytes/1024/1024); 

  return 0; 

} 
