/*
 引数 sm: numeric.js標準のCCS形式の疎行列
 引数 v:  密ベクトルの要素を並べた一次元配列
 戻り値:  smとvの積の密ベクトル
 numeric.jsでは疎行列とベクトルの積を計算するメソッドが用意されていないため自作
 密ベクトルをいったん疎行列形式に変換してccsDotを適用後，
 計算結果を密ベクトルに直して返す
*/
numeric.ccsDotMatVec = function ccsDotMatVec(sm, v) {
	var n = v.length;

	// 列インデックス [0, n, n, ... , n] (n+1要素)
	var colIdx = new Array(n + 1);
	colIdx[0] = 0;
	for(var i = 1; i < n + 1; i++)
		colIdx[i] = n;

	// 行インデックス [0, 1, 2, ... , n-1] (n要素)
	var rowIdx = new Array(n);
	for(var i = 0; i < n; i++)
		rowIdx[i] = i;

	// 疎行列をnumeric.jsの格納方式にならって作成
	var vSparseMat = [colIdx, rowIdx, v];

	var vecMat = numeric.ccsDot(sm, vSparseMat);

	// numeric.jsではCCS形式においてvalueが
	// 行インデックスの昇順に並ぶことが保障されないため
	// 非ゼロ値の並べ替えが必要
	var vec = new Array(n);
	for(var i=0; i<n; i++)
		vec[vecMat[1][i]] = vecMat[2][i];

	return vec;
}


/*
 引数 coo: COO形式の疎行列 coo[0]: 行Idx配列, coo[1]: 列Idx配列, coo[2]: 値配列
 戻り値: CCS形式の疎行列
 http://numericjs.com/workshop.php?link=4300ee0e171f8d0b7d316a2c1792b539b12a369ba4a3e3fef43bd114b0761aab
*/
numeric.cooToCcs = function cooToCcs(coo) {
	var B = numeric.sscatter(coo);
	var C = numeric.sgather(B);
	return numeric.ccsScatter(C);
}