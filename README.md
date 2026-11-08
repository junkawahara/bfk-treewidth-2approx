# Optimized 2-Approximation of Treewidth — C++ 実装

論文: M. Belbasi, M. Fürer, M. Kumar, *Optimized 2-Approximation of Treewidth*, arXiv:2411.16918v3
（PDF: `paper/2411.16918.pdf`）

グラフ G に対して、幅 ≤ 2k+1 の木分解を構成するか、tw(G) > k と結論する（Theorem 1）。
既定（自動）モードでは k を探索し、**幅 ≤ 2·tw(G)+1** の木分解と tw の下界を出力する。

## ビルドと実行

```
make                       # bin/tw2approx, bin/verify_td, bin/gen
./run_limited.sh ./bin/tw2approx graph.gr -o out.td     # 32GB / 10分 制限付き実行
./bin/verify_td graph.gr out.td                          # 木分解の検証
tests/battery.sh                                         # ランダムテスト（--check 付き）
```

`run_limited.sh` は `ulimit -v 32GB` と `timeout 600` をかけてコマンドを実行する。

### tw2approx のオプション

| オプション | 意味 |
|---|---|
| `-k K` | k を固定する。幅 ≤ 2K+1 の木分解を出力するか、`tw > K` を出力する（終了コード 1）。初期幅 ≤ 4K+3 が必要 |
| （`-k` なし） | 自動モード。Theorem 1 を繰り返し適用する（下記） |
| `--init mindeg\|random\|mixed` | 初期木分解を作る消去順序。min-degree（既定）／ランダム／確率 `--mix Q` でランダムな頂点を消去する min-degree（テスト用に、わざと悪い初期解を作る） |
| `--no-lb` | 退化度による下界 tw ≥ degeneracy を使わない（失敗経路のテスト用） |
| `--mem-gb G` | DP表のメモリ上限（既定 24GB）。超える場合は計算を始める前にエラー終了する |
| `--check` | Move・Split のたびに、grouped TD の不変条件と**全DP表を作り直した結果**を照合する（遅い・テスト用） |
| `-v` | ラウンドや k の試行を表示する |

入力は PACE `.gr` 形式、出力は PACE `.td` 形式（main ノードのみ出力する）。統計情報は stderr に出す。

`bin/gen` はテストグラフの生成器: `gen ktree N K P SEED`（部分 K-木で tw ≤ K）、`gen grid R C`、`gen gnm N M SEED`。

## 実装の対応表

| 論文 | 実装（`src/tw2approx.cpp`） |
|---|---|
| Def 6 grouped TD / Lemma 7 | `GTD::build`。消去順序から作る木分解はそのまま unique home（Lemma 5 の変換は不要）。intersection 子は bag ごとにまとめる |
| Def 14 DP表 `A_x` | `Node::E`。L^{|B|} 要素、bag 内の位置 p の頂点のラベルを L 進数の p 桁目に置く（0 = X, 1..3 = W_i） |
| 3.1 Initialization / Addition / Subtraction / Update | `legal_init`, `add_im`（main←inter）, `add_mi`（inter←main）, 符号 ±1 |
| 3.2 DFS instance | `dfs_find` |
| Alg 2 Move / Alg 3 Rotate | `move_to_root`, `rotate`（forget chain は `build_chain`） |
| Alg 4 Merge | `rotate` 内の merge と `attach_inter_*`（同じ bag の兄弟 inter ノードを合併する） |
| Def 10/11, Alg 5/6 Split | `split`, `process_pending`, `finish_split` |
| Remark 24（3分割と4分割） | ラウンドで \|B_r\| = w+1 ≥ 3k+4 なら L=3（C1,C2,S）、それ以外は L=4（C1,C2,C3,S）。切り替え時に全表を再計算する |

## 論文から補った点・修正した点

1. **DP表を「超過量」で持つ**: `E[P] = s(A[P]) − |X_P|`。こうすると Addition は
   main←inter が `E_m[P] += E_z[P|B_z]`、inter←main が `E_z[P] += min_l(E_y[P+v:l] + [l=X])` になる。
   Algorithm 5 の `s(A_xi[P']) = s(A_x[P]) − |Q|` は `E_xi[P'] = E_x[P]` になり、Merge も単純な加算で済む。
   どれも論文の式と等価である。

2. **Algorithm 5 の表の再利用式の修正（重要）**: 論文は、編集可能な子を引いた後の A_x を
   全コピー x_i に `s(A_x[P]) − |Q|` で流用している。しかしこの A_x には、**別のコピー x_j に移る非編集子 z
   の寄与**が残っている。B_z の X 部分に複数のラベルが割り当てられると、その寄与は 0 にならない。
   そのため x_i の表が過大評価になり、split を見落とす（誤って tw > k と判定する）おそれがある。
   本実装では、非編集子を最も多く受け取るコピー j\* だけでこの再利用式を使い、
   それ以外のコピーに移る非編集子は先に A_x から減算する。コピー j ≠ j\* の表は、
   Initialization の後に自分の非編集子を加算して作る。`--check` で、全表が作り直した結果と一致することを確認している。

3. **split の top-down 選択**: 編集可能な intersection ノード z の main 子 y で選べるのは、
   新しい頂点 v のラベルだけである。候補は {X} ∪ a(z) に限る（Def 11(b)）。
   a(z) にないラベルを a(z) 内のラベルにまとめても、分割は legal のままでサイズは増えない（Lemma 12 の議論）。
   そのため、制限しても Def 10 の「最小サイズ」は変わらない。サイズが同じなら X を優先する（|X| 最大）。
   根では、条件 `|W_i| + size < |B_r|`、`size ≤ k+1` を満たす分割のうち、サイズ最小、次に |X| 最大のものを選ぶ。

4. **非編集子の割り当て**: 非空の W が1つの非編集子は、そのラベルのコピーへ割り当てる。
   B_z がすべて X の非編集子は、「a(p(x)) の中の任意の1つ」（Def 11(a)）として j\* に割り当てる。

5. **コピーと Merge**: 非根ノード x のコピー x_i が新しい木で home を持つのは、x の home 頂点のラベルが i のときに限る。
   それ以外のコピーは祖父ノード側のコピーに含まれるので、作成せずに最初から merge する。
   その子や表への寄与は、最も近い実在する祖先コピーへ渡す。
   これは論文の「コピーを作ってから MERGE」と同じ結果になるが、無駄な表計算を省ける。

6. **DFS の再開**: 論文の再開規則（「編集されていない最初の pre-visited ノードから」）は曖昧である。
   そこで、split の後は根から DFS をやり直し、post-visited（closed）のノードとその部分木は飛ばす。
   closed な部分木は Move の経路にも編集領域にも入らないため、closed という印は split 後も有効である。

7. **初期木分解と k の探索**: 初期木分解は消去順序（min-degree）から作る。Corollary 2 の Bodlaender の縮約の代わりに、
   Theorem 1 を次のように繰り返し適用する。
   下界 lb（はじめは退化度）と現在の幅 cw を持ち、k = max(lb, ⌈(cw−3)/4⌉) で実行する。
   成功すると cw ≤ 2k+1 になる。失敗すると tw > k が分かるので lb = k+1 とする。cw ≤ 2·lb+1 になったら止める。
   このため初期幅が 4k+3 を超えていても 2-近似が保証される。

8. **連結成分ごとに処理する**: grouped TD は連結グラフに対して定義されている（根の bag S が空にならない）。
   そのため成分ごとに実行し、最後に木分解をつなぐ。

## 計算量とメモリの実際

DP表は全ノードに L^{|B|} 個の int32 を持つ。そのため、メモリはおおよそ n·L^{w+1}·4 バイト、
1回の表操作は L^{w+1} の時間がかかる（論文の O\*(81^k)n に対応する）。
実用上は初期幅がおよそ 16 以下（k ≤ 3〜4）が限界である。それを超える場合は、`--mem-gb` の上限を確認したうえで計算前にエラー終了する。

## テスト結果（このマシン、すべて 32GB / 10分 制限下）

- `tests/battery.sh`: ランダムな部分 K-木（K=1..3）、グリッド、G(n,m) の計 75 ケースを、`--check --init random --no-lb` で実行した。
  **失敗 0**、メモリガードで中断 10（ランダム初期化で初期幅が 19〜33 になったもの）。
  検証した内容: 出力が正しい木分解であること、幅 ≤ 2k+1、幅 ≤ 2·tw+1、下界 ≤ tw、全DP表が再計算結果と一致すること。
- `_GLIBCXX_DEBUG` と UBSan を有効にしたビルドで 25 ケース: 問題なし。
- 規模の大きい例: n=16000 の部分 3-木、初期幅 15（mixed 初期化）、`-k 3` で幅 7、2.1 秒、最大 RSS 354MB。
  n=2000〜8000 の部分 2-木、初期幅 9〜11、`-k 2` で幅 5、いずれも 0.05 秒未満。
