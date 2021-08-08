# DP-iso

An implement for [Efficient Subgraph Matching: Harmonizing DynamicProgramming, Adaptive Matching Order, and Failing Set Together (SIGMOD'19)](https://dl.acm.org/doi/abs/10.1145/3299869.3319880)

Deeply thankful to the code of [In-Memory Subgraph Matching: An In-depth Study(SIGMOD'20)](https://dl.acm.org/doi/abs/10.1145/3318464.3380581)



## Data Structure

* Adjacency List for Graph

  $[N_1(u_1),N_2(u1),...,N_1(u_2),N_2(u_2),...]$

where $N_i(u_j)$ for the ith neighbor of the jth vertex

by which both $G$ (data graph), $Q$ (query grpah) are stored .



* Candidates

  $C_1(u_1),C_2(u1),...,C_1(u2),C_2(u_2),...]$

where $C_i(u_j)$ for the ith candidate (in $G$) of the jth vertex in $Q$



## Filter Order

bfs tree of $Q$ + backward/forward edge definde by bfs order = DAG of $Q$

the start vertex of bfs is defined by $\text{argmin}_u\frac{|C(u)|}{degree(u)}$  , where $|C(u)|$ is the number of vertex u's candidates by LDF(Label and Degree Filter) 

```
1. StartVertex = argmin (|C(u)| / D(u))
2. BFSTree, BFSOrder = BFS(Q, StartVertex)
3. DAG = BuildDAG(BFSTree, BFSOrder)
```

## Filter

forward filter: filter order is defined by BFS order

backward filter: filter order is defined by reversed BFS order

DP-iso alternately use forward and backward filter to prune candidates



Take forward filter as an example,  DP-iso Filter try to prun all the $C(u)$ by 

```
Support = ForwardNeighbor(u)->C(u)
Count = |ForwardNeighbor(u)|
if Support < Count
	Prun(C(u))
```

## Index

* Edge Matrix $M$

  for each $M[u_1][u_2]$: $[Bitset(C_1),Bitset(C_2),...]$ 

  where $Bitset(C_i)$ is the Bitset of $C(u_2)$ related to $C_i(u_1)$ 

  Then given $u_1,u_2,M$, we can know if $C_i(u_1)$ is in the embedding, the valid candidates of $u_2$ is $M[u_1][u_2]$



## Matching Order

DP-iso use a DP of weight for mathing order, where the weight is an estimation of the maxinum embedding of paths.

$W(C(u_1)) = \min_{u_2} \sum_{C(u_2)}W(C(u_2))$ 

where $u_2$ is the backward neighbor / child of $u_1$, $C(u_2)$ is the valid candidates given $u_1$



Since the valid candidates change when enumarating, the matching order in DP-iso is adaptive, ad we use a priority queue to record the extendable vertices ($W(C(u)$) for first, $degree(u)$ for second) 

## Enumerating with Failing Set

DP-iso use a bitset for failing set, when enumerating/searching,

Searching at a leave node,

```
1. FailingSet = Empty, When success
2. FailingSet = Ancestor(u), When u has no candidates
3. FailingSet = Ancestor(u,v), When u,v confilct
```



For a non-leave node p,

```
1. FailingSet = Empty, if exist a Empty FailingSet in child node, when no failure exists
2. FailingSet = FailingSet(u), if child nodes u indicates the failure is not by p
3. FailingSet = Union(FailingSet(u)), for all chile nodes u, if child nodes u does not indicate the failure is not by p
```

