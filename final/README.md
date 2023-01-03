# Final exam

## Run
```
pip install biopython plotly dash pandas numpy gdown
python pipeline.py
```

## Q1

Count the number of variants in `twb2.csv` that the variant's ref is violate the reference in grch37.fa

(Ignore insertion)

```
{'ok': 734298, '-': 15671, 'violate': 2952}
```

## Q2

Count the number of variants `X->X` and plot two figures that decending with y value and ascending with x label.

(Ignore indel and strand)

```
ref alt    size
  A   C   28262
  A   G  104318
  A   T   10129
  C   A   36126
  C   G   18944
  C   T  148189
  G   A  150380
  G   C   19049
  G   T   33947
  T   A    9992
  T   C  106605
  T   G   26018
```


## Q3

Count the number of indel with it's length and plot it with insertion/deletion labels.

(Ignore multi-bases ref to mutli-bases alt, e.g. `AT->CC`)

```
     indel_len       type   size
0            1   deletion  20327
1            1  insertion  10600
2            2   deletion   8099
3            2  insertion   2004
4            3   deletion   2310
5            3  insertion    749
6            4   deletion   3643
7            4  insertion   1232
8            5   deletion   1233
9            5  insertion    313
10           6   deletion    336
11           6  insertion     81
12           7   deletion    522
13           7  insertion    121
14           8   deletion    434
15           8  insertion    113
16           9   deletion    172
17           9  insertion     51
18          10   deletion    278
19          10  insertion     46
20          11   deletion    276
21          11  insertion     39
22          12   deletion    137
23          12  insertion     27
24          13   deletion    183
25          13  insertion     37
26          14   deletion    155
27          14  insertion     43
28          15   deletion     73
29          15  insertion     16
30          16   deletion    104
31          16  insertion     20
32          17   deletion     72
33          17  insertion     21
34          18   deletion     42
35          18  insertion     16
36          19   deletion     70
37          19  insertion     23
38          20   deletion     57
39          20  insertion     14
40          21   deletion     45
41          21  insertion     10
42          22   deletion     38
43          22  insertion      7
44          23   deletion     37
45          23  insertion      8
46          24   deletion     17
47          24  insertion      1
48          25   deletion     26
49          25  insertion      2
50          26   deletion     27
51          26  insertion     10
52          27   deletion     23
53          27  insertion      3
54          28   deletion     23
55          28  insertion      5
56          29   deletion     20
57          29  insertion      5
58          30   deletion      7
59          30  insertion      1
60          31   deletion      3
61          31  insertion      6
62          32   deletion      8
63          33   deletion      5
64          34   deletion     10
65          34  insertion      4
66          35   deletion      9
67          35  insertion      1
68          36   deletion      7
69          37   deletion      5
70          37  insertion      1
71          38   deletion      4
72          39  insertion      1
73          40   deletion      2
74          40  insertion      4
75          41   deletion      4
76          41  insertion      1
77          42   deletion      3
78          43   deletion      2
79          43  insertion      1
80          44   deletion      5
81          45   deletion      2
82          46   deletion      3
83          46  insertion      2
84          47   deletion      5
85          47  insertion      1
86          48   deletion      1
87          49   deletion      4
88          50   deletion      2
89          50  insertion      1
90          51   deletion      1
91          52   deletion      4
92          52  insertion      2
93          53   deletion      3
94          53  insertion      1
95          54   deletion      3
96          54  insertion      1
97          55   deletion      1
98          55  insertion      1
99          56   deletion      3
100         57   deletion      2
101         57  insertion      1
102         59   deletion      1
103         59  insertion      2
104         60   deletion      1
105         61   deletion      1
106         62   deletion      3
107         63   deletion      1
108         64   deletion      2
109         64  insertion      3
110         65   deletion      3
111         66   deletion      2
112         67  insertion      2
113         68   deletion      2
114         70   deletion      2
115         71   deletion      6
116         72   deletion      7
117         73  insertion      1
118         74   deletion      1
119         75   deletion      1
120         76  insertion      1
121         78   deletion      4
122         80   deletion      1
123         81   deletion      2
124         82   deletion      2
125         84   deletion      2
126         86   deletion      1
127         87   deletion      2
128         88   deletion      2
129         89   deletion      2
130         92   deletion      1
131         93   deletion      1
132         94   deletion      5
133         96  insertion      2
134         97   deletion      1
135         98   deletion      1
136         99   deletion      1
137        103  insertion      1
138        105   deletion      2
139        106   deletion      3
140        108   deletion      2
141        109   deletion      2
142        110   deletion      1
143        110  insertion      1
144        113   deletion      1
145        114   deletion      2
146        115   deletion      2
147        121   deletion      1
148        122   deletion      3
149        123   deletion      2
150        126   deletion      4
151        127   deletion      1
152        128   deletion      1
153        129   deletion      1
154        131   deletion      1
155        133   deletion      2
156        134   deletion      1
157        136   deletion      1
158        137   deletion      1
159        138   deletion      2
160        143   deletion      1
161        144   deletion      2
162        146   deletion      2
163        147   deletion      3
164        149   deletion      2
165        152   deletion      3
166        154   deletion      1
167        155   deletion      2
168        155  insertion      1
169        159   deletion      1
170        160   deletion      1
171        161   deletion      2
172        162   deletion      1
173        163   deletion      1
174        164   deletion      1
175        165   deletion      2
176        166   deletion      3
177        170   deletion      1
178        171   deletion      3
179        174   deletion      1
180        175   deletion      1
181        176   deletion      1
182        177   deletion      1
183        191   deletion      2
184        192   deletion      1
185        193   deletion      1
186        194   deletion      2
187        199   deletion      1
188        200   deletion      2
189        200  insertion      1
190        202   deletion      2
191        204   deletion      1
192        205   deletion      1
193        212   deletion      1
194        215   deletion      1
195        236   deletion      2
196        243   deletion      1
197        246   deletion      1
198        248   deletion      1
199        249   deletion      2
200        253   deletion      2
201        259   deletion      1
202        264   deletion      1
203        268   deletion      1
204        269   deletion      1
205        270   deletion      1
206        270  insertion      1
207        280   deletion      1
208        281  insertion      2
209        283   deletion      2
210        284   deletion      2
211        300   deletion      1
212        309   deletion      1
213        319   deletion      1
214        329   deletion      1
215        337  insertion      1
216        338   deletion      2
217        339   deletion      1
218        345   deletion      2
219        363   deletion      1
220        378   deletion      1
221        388   deletion      1
222        392   deletion      1
223        409   deletion      1
224        412  insertion      1
225        421   deletion      1
226        436   deletion      1
227        446   deletion      2
228        450   deletion      1
229        454  insertion      2
230        507   deletion      1
231        509   deletion      1
232        510   deletion      2
233        543   deletion      1
234        559   deletion      1
235        606   deletion      1
236        645   deletion      1
237        649   deletion      1
238        657   deletion      1
239        679   deletion      1
240        709   deletion      2
241        748   deletion      1
242        750   deletion      1
243        753   deletion      2
244        768   deletion      1
245        773   deletion      1
246        780   deletion      1
247        806   deletion      1
248        808   deletion      1
249        898   deletion      1
250        966   deletion      3
251        984   deletion      1
252        988   deletion      2
253       1008   deletion      2
254       1037   deletion      3
255       1051   deletion      2
256       1058  insertion      2
257       1066   deletion      1
258       1072   deletion      1
259       1084   deletion      1
260       1102   deletion      2
261       1123   deletion      1
262       1137   deletion      2
263       1158   deletion      1
264       1378   deletion      1
265       1403   deletion      2
266       1458   deletion      2
267       1563   deletion      2
268       1694   deletion      1
269       1759   deletion      1
270       1790   deletion      2
271       1849   deletion      1
272       1871   deletion      1
273       1879   deletion      2
274       1885   deletion      2
275       2031   deletion      1
276       2046   deletion      2
277       2051   deletion      2
278       2118   deletion      1
279       2179   deletion      1
280       2269   deletion      2
281       2306   deletion      1
282       2348   deletion      1
283       2365   deletion      1
284       2406   deletion      1
285       2409   deletion      1
286       2410   deletion      1
287       2411   deletion      1
288       2631   deletion      1
289       2680   deletion      2
```
