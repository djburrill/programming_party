/*
 OpenMP test

 Calculate value of pi through parallel integration
*/

#include <omp.h>
#include <iostream>

int main(){
    // Set number of threads
    omp_set_dynamic(0);
    omp_set_num_threads(800);

    // Variables
    int numThreads = omp_get_max_threads();
    std::cout << "Number of Threads: " << numThreads << std::endl;
    int numSteps = 4;

    double sumVal = 0;
    double sumArray [numThreads*numSteps];
    double threadArray [numThreads*numSteps];
    double binSize = 1.0/numThreads;

    // Use OMP to parallelize
    #pragma omp parallel
    {
        int threadNum = omp_get_thread_num();
        double xVal = 0.0;

        // Define endpoints and step size
        double xStart = (threadNum)*binSize;
        double xEnd = (threadNum+1)*binSize;
        double stepSize = (xEnd-xStart)/numSteps;

        for (int i=0; i<numSteps; i++){
            xVal = ((xStart+i*stepSize)+(xStart+(i+1)*stepSize))/2.0;
            sumArray[threadNum*numSteps+i] = (4.0/(1+(xVal*xVal)))*stepSize;
            threadArray[threadNum*numSteps+i] = threadNum;
        }
    }

    // Sum up values of array
    for (int i=0; i<numThreads*numSteps; i++){
        sumVal += sumArray[i];
    }

    // Print out calculated value
    std::cout << "Pi: " << sumVal << std::endl;
}
