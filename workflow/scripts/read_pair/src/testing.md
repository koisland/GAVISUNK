```python
kmer_set = set()
for k in records:
    reads_dis = k[-2]
    ctg_dis = k[-1]

    if  abs(ctg_dis - reads_dis) / reads_dis < dis_threshold:
        kmer_set.add(k[0])
        kmer_set.add(k[1])
```

```python
valid_pairs = 0
for k in records:
    reads_dis = k[-2]
    ctg_dis = k[-1]

    if  abs(ctg_dis - reads_dis) / reads_dis < dis_threshold:
        valid_pairs += 1
```