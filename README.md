# Artificial Inteligence Project Report
Author: Harkeerat Singh Sawhney

## Summery of the project:
I have chosen to implement the Traveling Sales Man Problem with Ant Colony Optimization. I also used Nearest Neighbor as my initial starting result, and then ran Ant Colony on the results which I got from that. The main aim for this was so that I can get myself into a good starting position. Alongside that I implemented 2-opt as my local search algorithm. After I find my solution I run 2-opt on it to make the solution better. I tried running 2-opt inside Ant Colony Itself, however the implementation of mine was not giving my the correct results which I wanted. I faced this bug towared the last minute, which led me to spend a tough all night to fix this and get the correct results. The results of the project are attached in the excel sheet, which was given to us.

For the bigger problems, I had to use 2-opt right after the nearest Neightbor, as othweise my results were very far of. With the improved of Nearest neighbor with 2-opt it was able to give me a good starting point to run my Ant Colony over.

I gave more time for 2-opt to run as it improved my results drasticly. Two-opt gave me best results when I ran it after short problems, and before and after long problems.

Overall my average gap for all the problems are **10%**

## Problems I faced:
The biggest problem I faced was the bug in my 2-opt implementation. I was not able to get the correct results, which was very frustrating. I spent a lot of time trying to fix this, and I was able to fix it in the end. I also faced a lot of problems with the Ant Colony, as I was not able to get the correct results. I was able to fix this by changing my parameters and also by changing the way I was calculating the pheromone. Along with that in Ant Colony I had the problem with my solution not being valid. It was difficult to debugg, since unlike python, I couldnt not graph it.

## Possible Improvements
My results can deffinatly be improved. I feel I can fix the Two-opt I can achieve much better results, since I can run this inside the Ant-colony itself. Also if I ran this multiple times to get a better starting point, with my seed and also to improve on my parameters, it could have improved the score a lot.

## How to run the project:
In order to run the project the following commmands need to be run:
```
g++ -Ofast ./TSP_ACP.cpp -o tsp_acp && ./tsp_acp
```
