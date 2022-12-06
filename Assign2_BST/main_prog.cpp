/*
  CS6013: Advanced Data Structures and Algorithms
  Programming Assignment II

Submitted By,
    Shrusti
    CS22MTECH11017
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <iterator>
#include <cfloat>
#include <cstring>

using namespace std;

float dp[101][101];
int root[101][101]; // table storing the root of the tree which gives the minimum cost

// find the sum of probabilities
float sumOfProb(float prob[], int i, int j)
{
    float sum = 0.0;
    for (int k = i; k < j; k++)
    {
        sum += prob[k];
    }
    return sum;
}

// function to find optimal BST using Dynamic Programming
float optimalBST(float prob[], int n)
{
    // if single key, then probOfOccurence is the cost
    for (int i = 0; i < n; i++)
    {
        dp[i][i + 1] = prob[i];
    }

    // if number of nodes is > 1 i.e; 2,3, 4...
    for (int len = 2; len <= n + 1; len++)
    {
        // row of dp[][]
        for (int i = 0; i <= n - len + 1; i++)
        {
            // j denotes column number
            int j = i + len - 1;

            dp[i][j] = FLT_MAX;

            // sum of probabilties from i to j
            float weight = sumOfProb(prob, i, j);

            // trying to make all nodes as root of BST
            for (int k = i + 1; k <= j; k++)
            {
                // cost when k becomes the root
                float costOfTree;

                // left tree: nodes from i to (k-1) are on left of the tree
                float leftTree = dp[i][k - 1];

                // right tree : nodes from k to j are on right of the tree
                float rightTree = dp[k][j];

                // calculating overall cost
                costOfTree = leftTree + rightTree + weight;

                // if minimum cost is found, update the dp table and also store the root of the tree
                if (costOfTree < dp[i][j])
                {
                    dp[i][j] = costOfTree;
                    root[i][j] = k;
                }
            }
        }
    }
    return dp[0][n];
}

// function to print preorder traversal of the Optimal BST
void preOrderTraversal(int i, int j, string input[])
{
    if (i == j)
        return;
    // obtaining the index of the root
    int key = root[i][j];
    // printing the root
    cout << input[key - 1] << " ";

    // recursive calls for left and right subtree, until you reach leaf nodes of the tree
    preOrderTraversal(i, key - 1, input);
    preOrderTraversal(key, j, input);
}

int main()
{
    // variable declaration
    int n;

    cout << "How many strings do you want to insert in the BST ?" << endl;
    cin >> n;
    string input[n];
    float prob[n];

    cout << "Enter " << n << " strings in sorted dictionary order along with their probabilities :" << endl;
    // for taking input
    for (int i = 0; i < n; i++)
    {
        cin >> input[i];
        cin >> prob[i];
    }
    // check if input is sorted or not
    bool isInputSorted = is_sorted(input, input + n);
    if (!isInputSorted)
    {
        cout << "The strings entered are not in sorted order." << endl;
        return 0;
    }

    // check if duplicate probabililties exist using Set
    unordered_set<float> distinctProb;
    for (int i = 0; i < n; i++)
    {
        distinctProb.insert(prob[i]);
    }
    if (n != distinctProb.size())
    {
        cout << "The probabilities are not distinct." << endl;
        return 0;
    }

    // check if probabilities sum to one or not
    float sum = sumOfProb(prob, 0, n);
    if (sum < 1 || sum > 1)
    {
        cout << "The probabilities dont add up to 1." << endl;
        return 0;
    }

    // Find the BST that provides the minimum expected total access time.
    float accessTime = optimalBST(prob, n);
    cout << endl<<"The minimum expected total access time is : " << accessTime << endl;

    // print a preorder traversal of this BST.
    cout << "Preorder traversal of the BST that provides minimum expected total access time is: " << endl;
    preOrderTraversal(0, n, input);
    cout << endl;
}