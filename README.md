# HOQRI: Higher-order QR Iteration for Scalable Tucker Decomposition
Yuchen Sun, Kejun Huang

## Citing

## Requirements
C++

## Input file format
The data file lists all nonzeros elements. The first three colums are the indexes and the last colume is the corresponding value.
Take a 5\*5\*5 tensor for example, the number of nonzero elements is 6. The tensor example exists in the root directory, named 'mytensor.txt':
>
<dl><pre>
1 1 2 1.5
1 2 2 2.5
2 1 1 3.7
1 2 3 0.5
2 1 2 4.1
5 4 4 1.1
</pre></dl>

The Facebook, 1M-MovieLen, Last and Delicious datasets are in the 'dataset.zip'. 10M-MovieLen is not uploaded due to the large size.

As shown in our paper, the dimensions of these datasets are
| Dataset  | Size | nnz |
| ------------- | ------------- | ------------- |
| Facebook  | 46952, 46951, 1592  | 738079 |
| MovieLen | 610, 49961, 8215 | 84159 |
| Last | 2100, 18744, 12647 | 186479 |
| Delicious | 108035, 107253, 52955 | 437593 |
| 10M-MoiveLen | 162541, 49994, 9083 | 20503478 |

## How to run
Go to the path of downloaded project

#### Compile
Linux and Windows system: `g++ *.cpp -o HOQRI`

MacOS: `clang++ -std=c++11 *.cpp -o HOQRI`

#### Run 
`<install_path>/HOQRI <TensorDimension> <Rank> <DataFilePath>`, TensorDimension and Rank are of size 3, splitted by comma.

For example, `<install_path>/HOQRI 5,5,5 2,3,4 '<install_path>/mytensor.txt'`
