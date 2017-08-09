table.cpp & table.h are the files to implement pagerank algorithm. <br />
final.cpp - This file is responsible for the complete flow of this project. <br />
<br />
To generate the graph instances for different values of edge weights(probability), uncomment the line no. 291 in final.cpp. <br />
final.cpp uses the optimised algorithm. If there is a need to run the greedy version of the algorithm, then change the word 'celf' at line 295 to 'greedy'. <br />
<br />
Finally run on the terminal using the following command:<br />
> g++  -std=c++11 final.cpp <br />
> ./a.out    <br />
