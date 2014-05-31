// JavaScript Document
/// <reference path="http://ajax.googleapis.com/ajax/libs/jquery/1/jquery.min.js" />
/// <reference path="numeric-1.2.6.min.js" />
/// <reference path="outline.js" />
/// <reference path="delaunay.js" />
/// <reference path="fem.js" />



////////////////////////////////////////////////////////////
// FEMクラスの定義
////////////////////////////////////////////////////////////
	
function FEM(initpos, tri){
	this.pos = numeric.clone(initpos);      // 節点現在位置
	this.initpos = numeric.clone(initpos);  // 節点初期位置

	this.pos.push([0,0]);	// ゴミを追加、こいつを固定点にすることでなぜか解が安定する
	this.initpos.push([0,0]);
	this.posNum = this.pos.length;

	this.tri = numeric.clone(tri); // 三角形要素の節点リスト
	this.nodeToTri = [];		// ノードに接続している三角形要素リスト
	this.surEdge = [];			// 表面エッジリスト
	this.surToTri = [];		// 表面エッジ-対応する三角形要素リスト
	this.triToSur = [];		// 三角形要素-対応する表面エッジリスト
	this.surNode = [];		// 表面頂点リスト
	this.sndToSur = [];	// 表面頂点に接続している表面エッジリスト
	this.sndToUnDupTri = [];	// 表面頂点に接続している三角形のうち表面エッジを有していないもののリスト
	this.colNdFlag = numeric.rep([this.posNum], 0);	// for visualization
	this.colTriFlag = numeric.rep([this.triNum],0);
	this.normal = [];	// 表面エッジの法線ベクトル
	this.ndNormal = [];	// 頂点法線ベクトル


	this.makeSurface();


	this.Re = [];           // 要素回転マトリクス
	this.Be = [];           // 要素 ひずみマトリクス
	this.De = [];           // 要素 ひずみ-応力変換マトリクス
	this.Ke = [];           // 要素剛性マトリクス    
	var young = 100;       // ヤング率 [Pa]
	var poisson = 0.3;      // ポアソン比
	var density = 0.001;    // 密度 [kg/mm3]
	var thickness = 1;  // 物体の厚さ [mm]
	this.mass = [];     // 節点の等価質量
	this.alpha = 0.02;  // Mに作用するレイリー減衰のパラメータ
	this.beta = 0.01;    // Kに作用するレイリー減衰のパラメータ
	this.gravity = 100;

	this.penalty = 100;	// ペナルティ法による自己接触の係数

    // 要素剛性マトリクス作成
	this.makeMatrixKe(young, poisson, density, thickness);

	this.K = [];
	this.makeMatrixK();
	this.Vel = numeric.linspace(0,0,2*this.posNum);
	this.th = numeric.linspace(0,0,this.tri.length);
	this.foffset = numeric.linspace(0,0,2*this.posNum);
	this.flist = [];
	this.dlist = [];
	this.f = [];
	this.ff = [];
	this.ud = [];
	this.uf = [];		
	this.vf = [];
	this.vd = [];
	this.vdbuf = [];
	this.Kff = [];
	this.Kfd = [];
	this.Kdf = [];
	this.Kdd = [];
	this.Mf = [];
	this.Md = [];
	this.xd = [];
	this.xf = [];
	this.fof = [];
	this.fod = [];


	this.maxPStress = [];
	this.fixNode = []; // 固定ノードリスト: パブリックアクセスで外部から設定する
    // 破壊に関するメンバ変数
    this.removedFlag = new Array(this.tri.length);
    for(var i=0; i<this.tri.length; i++)
        this.removedFlag[i] = false;
    this.thrPStress = 5*young;

	// マウスでつまむためのメンバ変数
	this.holdNode = [];
	this.mousePosClick = [];
	this.uClick = []; // setBoundaryのためのメンバ, クリック時のUベクトル
	this.gripRad = 50; // setBoudaryにおける周辺拘束領域の半径
}


FEM.prototype.makeSurface = function(){
	// nodeToTriの作成
	this.nodeToTri = new Array(this.posNum);
	for(var i=0; i<this.posNum; i++)
		this.nodeToTri[i] = [];
	for(var i = 0; i < this.tri.length; i++) 
		for(var vert = 0; vert < 3; vert++) 
			this.nodeToTri[this.tri[i][vert]].push(i);

	// surEdge, surToTri, triToSur, sndToSurの作成
	// 四面体についてのループを回し、
	// 四面体の各エッジが現在着目している四面体以外の四面体に共有されているかどうかを調べる
	// (エッジの頂点番号からnodeToTetを参照すれば判定できる)
	// 共有されていなければそれは表面エッジであるとみなす
	var buf = [[0,1,2],[1,2,0],[2,0,1]];
	var n1,n2,n3;	// bufに対応する頂点番号
	var nt1, nt2; // nodeToTetの一次格納変数
	var shareFlag;
	var surCount = 0;
	var v1,v2;	
	this.triToSur = new Array(this.tri.length);
	for(var i = 0; i < this.tri.length; i++) {
		this.triToSur[i] = [];
	}
	for(var i = 0; i < this.tri.length; i++) {
		for(var edg = 0; edg < 3; edg++) {
			shareFlag = false;
			n1 = this.tri[i][buf[edg][0]];
			n2 = this.tri[i][buf[edg][1]];
			nt1 = this.nodeToTri[n1];
			nt2 = this.nodeToTri[n2];
			for(var j = 0; j < nt1.length; j++) {
				for(var k = 0; k < nt2.length; k++) {
					if(shareFlag)break;
					if(nt1[j] === nt2[k] && nt1[j] !== i) {
						shareFlag = true;
					}
				}
			}
			if(!shareFlag) {
				// surEdgeに格納する頂点番号の順番が反時計回りになるようにする
				n3 = this.tri[i][buf[edg][2]];
				v1 = numeric.sub(this.initpos[n1],this.initpos[n3]);
				v2 = numeric.sub(this.initpos[n2],this.initpos[n3]);
				if(v1[0]*v2[1]-v1[1]*v2[0]>0)
					this.surEdge.push([this.tri[i][buf[edg][0]], this.tri[i][buf[edg][1]]]);
				else
					this.surEdge.push([this.tri[i][buf[edg][1]], this.tri[i][buf[edg][0]]]);
				this.surToTri.push(i);
				this.triToSur[i].push(surCount);
				++surCount;
			}
		}
	}

	// surNode, sndToSurの作成
	this.surNode = [];
	var nd, dupFlag;
	for(var edg = 0; edg < this.surEdge.length; edg++) {
		for(var i = 0; i < 2; i++) {
			dupFlag = false;
			nd = this.surEdge[edg][i];
			for(var j = 0; j < this.surNode.length; j++) {
				if(nd === this.surNode[j]) {
					dupFlag = true;
					this.sndToSur[j].push(edg);
					break;
				}
			}
			if(!dupFlag){
				this.surNode.push(nd);
				this.sndToSur.push([edg]);
			}
		}
	}

	// sndToUnDupTriの作成
	this.sndToUnDupTri = new Array(this.surNode.length);
	for(var i=0; i<this.surNode.length; ++i)
		this.sndToUnDupTri[i] = [];
	var nd;
	var dupFlag;
	for(var snd=0; snd<this.surNode.length; ++snd){
		nd = this.surNode[snd];
		for(var ntri=0; ntri<this.nodeToTri[nd].length; ++ntri){
			dupFlag = false;
			// 着目する三角形がsurToTriの中のすべての三角形と重複がなければ
			// 内側三角形リストに追加する
			for(var ond = 0; ond < this.surNode.length; ++ond) {
				for(var edg=0; edg<this.sndToSur[ond].length; ++edg){
					if(this.nodeToTri[nd][ntri] == this.surToTri[this.sndToSur[ond][edg]]){
						dupFlag = true;
						break;
					}
				}
				if(dupFlag)break;
			}
			if(!dupFlag)
				this.sndToUnDupTri[snd].push(this.nodeToTri[nd][ntri]);
		}
	}
}

FEM.prototype.makeMatrixKe = function(young, poisson, density, thickness){
	// Bマトリクスを作成
	var TriNum = this.tri.length;
	this.Be = new Array(TriNum);
	for(var i=0; i<TriNum; i++){
		this.Be[i] = makeMatrixB(this.pos[this.tri[i][0]], this.pos[this.tri[i][1]], this.pos[this.tri[i][2]]);
	}
	// Dマトリクスを作成
	this.De = new Array(TriNum);
	for (var i = 0; i < TriNum; i++) {
		this.De[i] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
		var tmp = young / (1.0 - poisson * poisson);
		this.De[i][0][0] = tmp;
		this.De[i][0][1] = poisson * tmp;
		this.De[i][1][0] = poisson * tmp;
		this.De[i][1][1] = tmp;
		this.De[i][2][2] = 0.5 * (1 - poisson) * tmp;
	}
	// Keマトリクスを作成
	var KeTmp;
	var Bt;
	var posMat;
	var area;
	this.mass = numeric.linspace(0,0,this.pos.length);
	for(var i=0; i<TriNum; i++){
		Bt = numeric.transpose(this.Be[i]);
		KeTmp = numeric.dot(this.De[i],this.Be[i]);
		KeTmp = numeric.dot(Bt,KeTmp);
		posMat =  [
		[1,this.pos[this.tri[i][0]][0],this.pos[this.tri[i][0]][1]], 
		[1,this.pos[this.tri[i][1]][0],this.pos[this.tri[i][1]][1]], 
		[1,this.pos[this.tri[i][2]][0],this.pos[this.tri[i][2]][1]] ];
		area = 0.5 * numeric.det(posMat);
		area = Math.abs(area);
		KeTmp = numeric.mul(KeTmp,area*thickness);
		this.Ke.push(KeTmp);
		
		for(var j=0; j<3; j++)
			this.mass[this.tri[i][j]] += area * density * thickness * 0.333333;			
	}
    // stress, maxPStressの初期化
	this.stress = new Array(TriNum);
	for (var i = 0; i < TriNum; i++) {
	    this.stress[i] = [0, 0, 0];
	}
	this.maxPStress = numeric.linspace(0, 0, TriNum);

    // Reの初期化
	this.Re = new Array(TriNum);
	for(var i=0; i<TriNum; i++){
	    this.Re[i] = numeric.identity(2);
	}
}

// 三角形定ひずみ要素の全体剛性マトリクスを作成する関数
FEM.prototype.makeMatrixK = function(){
	this.K = numeric.rep([2*this.pos.length,2*this.pos.length],0);
	for(var i=0; i<this.tri.length; i++){
		for(var j=0; j<3; j++){
			for(var k=0; k<3; k++){
				for(var l=0; l<2; l++){
					for(var m=0; m<2; m++){
						this.K[2*this.tri[i][j]+l][2*this.tri[i][k]+m] += this.Ke[i][2*j+l][2*k+m];
					}
				}
			}
		}
	}
}


// クリック時の処理
FEM.prototype.selectHoldNodes = function(mousePos){	

	this.mousePosClick = new Array(mousePos.length);
	for(var i=0; i<mousePos.length; i++){
		this.mousePosClick[i] = new Array(2);
		this.mousePosClick[i][0] = mousePos[i][0];
		this.mousePosClick[i][1] = mousePos[i][1];
	}
		
	
	
	this.uClick = new Array(mousePos.length);
	for(var i=0; i<mousePos.length; i++){
		this.uClick[i] = numeric.linspace(0,0,2*this.pos.length);
	}
	
	this.holdNode = new Array(mousePos.length);
	for(var i=0; i<mousePos.length; i++){
		this.holdNode[i] = [];				
	}
	
	var dif;
	var dist;
	for(var cl=0; cl<mousePos.length; cl++){
		for(var i=0; i<this.pos.length; i++){
			dif = numeric.sub(mousePos[cl],this.pos[i]);
			dist = numeric.norm2(dif);
			if(this.gripRad > dist){
				this.uClick[cl][2*i] = this.pos[i][0]-this.initpos[i][0];
				this.uClick[cl][2*i+1] = this.pos[i][1]-this.initpos[i][1];
				this.holdNode[cl].push(i);
			}
		}
	}
}


// 境界条件の設定
FEM.prototype.setBoundary = function(clickState, mousePos, gravityFlag, selfCollisionFlag){
	
	if(mousePos.length != this.holdNode.length)
		this.selectHoldNodes(mousePos);
	
	this.dlist = [];
	this.flist = [];
	this.ud = [];
	this.ff = [];
	
	var nodeToDF = numeric.linspace(0,0,this.posNum);
	var u = numeric.linspace(0,0,2*this.posNum);
	this.f = numeric.linspace(0,0,2*this.posNum);

	// 固定ノードの境界条件
	for(var i=0; i<this.fixNode.length; i++) {
		var nd=this.fixNode[i];
		u[2*nd] = this.pos[nd][0]-this.initpos[nd][0];
		u[2*nd+1]=this.pos[nd][1]-this.initpos[nd][1];
		nodeToDF[nd]="d";
	}
	
	// 上面のノードを固定
	for(var i=0; i<this.pos.length; i++){
		if(i==this.pos.length-1){
			u[2*i] = 0;
			u[2*i+1] = 0;
			nodeToDF[i] = "d";
		}else if(nodeToDF[i]!="d"){
			this.f[2*i] = 0;
			if(gravityFlag)
				this.f[2*i+1] = this.gravity * this.mass[i];
			else
				this.f[2*i+1] = 0;
			
			nodeToDF[i] = "f";
		}
	}
	
	// クリックノードの境界条件
	if(clickState == "Down"){
		for(var cl=0; cl<mousePos.length; cl++){
			for(var i=0; i<this.holdNode[cl].length; i++){
				var nd = this.holdNode[cl][i];
				if(nodeToDF[nd]=="d")continue;
				u[2*nd]   = this.uClick[cl][2*nd]+mousePos[cl][0]-this.mousePosClick[cl][0];
				u[2*nd+1] = this.uClick[cl][2*nd+1]+mousePos[cl][1]-this.mousePosClick[cl][1];
				nodeToDF[nd] = "d";
			}
		}
	}

	// エッジ法線ベクトルの作成
	this.normal = numeric.rep([this.surEdge.length, 2], 0);
	var pe0, pe1;	// エッジの頂点位置ベクトル
	var veclen;
	var normalTmp;
	for(var sur = 0; sur < this.surEdge.length; ++sur) {
		pe0 = this.pos[this.surEdge[sur][0]];
		pe1 = this.pos[this.surEdge[sur][1]];
		normalTmp = [pe1[1] - pe0[1], -(pe1[0] - pe0[0])];
		veclen = numeric.norm2(normalTmp);
		if(veclen===0) continue; // 法線ベクトルの長さがゼロになる場合は前回までの値を採用する
		this.normal[sur] = numeric.div(normalTmp, veclen);
	}

	// 頂点法線ベクトルの作成
	this.ndNormal = numeric.rep([this.surNode.length,2],0);
	var ndNormalTmp, ndNmNorm;
	for(var snd = 0; snd < this.surNode.length; ++snd) {
		ndNormalTmp = [0,0];
		for(var sedg = 0; sedg < this.sndToSur[snd].length; ++sedg) 
			ndNormalTmp = numeric.add(ndNormalTmp, this.normal[this.sndToSur[snd][sedg]]);
		ndNmNorm = numeric.norm2(ndNormalTmp);
		if(ndNmNorm===0)continue; // 法線ベクトルの長さがゼロになる場合は前回までの値を採用する
		this.ndNormal[snd] = numeric.div(ndNormalTmp, ndNmNorm);
	}

	// 三角形の内外判定によるペナルティ法
	if(selfCollisionFlag) {
		var nd;		// 着目するノード番号
		var tr;		// 着目する三角形番号
		var p0,p1,p2;	// 三角形の頂点位置ベクトル
		var pe0;	// エッジの頂点位置ベクトル
		var q;			// 着目するノードの位置ベクトル
		var v0,v1,v2;	// 辺ベクトル
		var q0,q1,q2;	// i番頂点からqまでの相対位置ベクトル
		var normal;		// 表面エッジの法線ベクトル
		var qd;			// エッジ0から着目するノードへの相対位置ベクトル
		var d;			// くいこみ量
		// for visualization
		this.colNdFlag = numeric.rep([this.posNum],0);
		this.colTriFlag = numeric.rep([this.tri.length],0);
		for(var sur = 0; sur < this.surEdge.length; sur++) {
			tr = this.surToTri[sur];
			p0 = this.pos[this.tri[tr][0]];
			p1 = this.pos[this.tri[tr][1]];
			p2 = this.pos[this.tri[tr][2]];
			pe0 = this.pos[this.surEdge[sur][0]];
			// 辺ベクトルの作成
			v0 = numeric.sub(p1,p0);
			v20 = numeric.sub(p2,p0);
			v1 = numeric.sub(p2,p1);
			v2 = numeric.sub(p0,p2);
			normal = this.normal[sur];
			for(var i = 0; i < this.surNode.length; i++) {
				// 着目するノードが三角形要素の頂点なら無視
				nd = this.surNode[i];
				if(this.tri[tr][0]===nd)continue;
				if(this.tri[tr][1]===nd)continue;
				if(this.tri[tr][2]===nd)continue;
				// ノードの三角形内外判定
				q = this.pos[nd];
				q0 = numeric.sub(q,p0);
				q1 = numeric.sub(q,p1);
				q2 = numeric.sub(q,p2);
				// 三角形要素の格納順番が時計回りか反時計回りかによって
				// 三角形に対する頂点の内外判定の符号を変える
				if(v0[0] * v20[1] - v0[1] * v20[0] > 0) {
					if(v0[0]*q0[1]-v0[1]*q0[0]<0 
						|| v1[0]*q1[1]-v1[1]*q1[0]<0 
						|| v2[0]*q2[1]-v2[1]*q2[0]<0 )
						continue;
				} else {
					if(v0[0]*q0[1]-v0[1]*q0[0]>0 
						|| v1[0]*q1[1]-v1[1]*q1[0]>0 
						|| v2[0]*q2[1]-v2[1]*q2[0]>0 )
						continue;
				}
				// くいこみ量の計算
				qd = numeric.sub(q, pe0);
				d = numeric.dot(normal,qd);
				if(isNaN(normal[0])) alert("nm0");
				if(isNaN(normal[1])) alert("nm1");
				if(isNaN(d)) alert("nan1");
				// 外力の設定
				for(var vt = 0; vt < 2; vt++) {
					this.f[2*this.surEdge[sur][vt]+0] += this.penalty*d*normal[0]*0.5;
					this.f[2*this.surEdge[sur][vt]+1] += this.penalty*d*normal[1]*0.5;
				}
				// 頂点への反作用
				this.f[2*nd+0] += -this.penalty*d*normal[0];
				this.f[2*nd+1] += -this.penalty*d*normal[1];
				// for visualization
				this.colNdFlag[nd]=1;
				this.colTriFlag[tr]=1;
			}
		}

		// 表面内側の三角形との接触判定
		for(var snd=0; snd<this.surNode.length; ++snd){
			for(var str=0; str<this.sndToUnDupTri[snd].length; ++str){
				tr = this.sndToUnDupTri[snd][str];
				p0 = this.pos[this.tri[tr][0]];
				p1 = this.pos[this.tri[tr][1]];
				p2 = this.pos[this.tri[tr][2]];
				// 辺ベクトルの作成
				v0 = numeric.sub(p1,p0);
				v20 = numeric.sub(p2,p0);
				v1 = numeric.sub(p2,p1);
				v2 = numeric.sub(p0,p2);
				normal = this.ndNormal[snd];
				for(var i = 0; i < this.surNode.length; i++) {
					// 着目するノードが三角形要素の頂点なら無視
					nd = this.surNode[i];
					if(this.tri[tr][0]===nd)continue;
					if(this.tri[tr][1]===nd)continue;
					if(this.tri[tr][2]===nd)continue;
					// ノードの三角形内外判定
					q = this.pos[nd];
					q0 = numeric.sub(q,p0);
					q1 = numeric.sub(q,p1);
					q2 = numeric.sub(q,p2);
					// 三角形要素の格納順番が時計回りか反時計回りかによって
					// 三角形に対する頂点の内外判定の符号を変える
					if(v0[0] * v20[1] - v0[1] * v20[0] > 0) {
						if(v0[0]*q0[1]-v0[1]*q0[0]<0 
							|| v1[0]*q1[1]-v1[1]*q1[0]<0 
							|| v2[0]*q2[1]-v2[1]*q2[0]<0 )
							continue;
					} else {
						if(v0[0]*q0[1]-v0[1]*q0[0]>0 
							|| v1[0]*q1[1]-v1[1]*q1[0]>0 
							|| v2[0]*q2[1]-v2[1]*q2[0]>0 )
							continue;
					}
					// くいこみ量の計算
					qd = numeric.sub(q, this.pos[this.surNode[snd]]);
					d = numeric.dot(normal,qd);
					// 外力の設定
					// 食い込んだ頂点への反発力
					this.f[2*this.surNode[snd]+0] += this.penalty*d*normal[0];
					this.f[2*this.surNode[snd]+1] += this.penalty*d*normal[1];
					// 頂点への反作用
					this.f[2*nd+0] += -this.penalty*d*normal[0];
					this.f[2*nd+1] += -this.penalty*d*normal[1];
					// for visualization
					this.colNdFlag[nd]=2;
					this.colTriFlag[tr]=2;
				}				
			}
		}
	}

	for(var i=0; i<this.pos.length; i++){
		if(nodeToDF[i] == "d"){
			this.dlist.push(i);
			this.ud.push(u[2*i]);
			this.ud.push(u[2*i+1]);
		}else{
			this.flist.push(i);
			this.ff.push(this.f[2*i]);
			this.ff.push(this.f[2*i+1]);
		}
	}
}



// Stiffness Warpe 法のためのKマトリクスを作成する
// modeはMovingAverageかNormal
FEM.prototype.makeMatrixKSW = function(){
	var TriNum = this.tri.length;
	this.K = numeric.rep([2*this.pos.length,2*this.pos.length],0);
	var RK = numeric.rep([2*this.pos.length,2*this.pos.length],0);
	
	for(var i=0; i<TriNum; i++){
        // 削除された要素は計算しない
        if(this.removedFlag[i]) continue;
		
		// 要素の回転角度取得
		var q = numeric.rep([3,2],0);
		for(var j=0; j<2; j++)
			q[j] = numeric.sub(this.initpos[this.tri[i][j+1]],this.initpos[this.tri[i][0]]);
		
		var p = numeric.rep([3,2],0);
		for(var j=0; j<2; j++)
			p[j] = numeric.sub(this.pos[this.tri[i][j+1]], this.pos[this.tri[i][0]]);
		
		var A = 0;
		for(var j=0; j<2; j++)
			A += q[j][1]*p[j][0] - q[j][0]*p[j][1];
		
		var B = 0;
		for(var j=0; j<2; j++)
			B += q[j][0]*p[j][0] + q[j][1]*p[j][1];

		// 回転角度の計算
		var th_now = -Math.PI/2.0+Math.atan2(B,A);
		
		// 安定化のために回転場の移動平均
		// atan2の不連続性を考慮して前回との開き角度が180度を超えた場合は
		// 鋭角側の差をとることにする
		/*
		if(Math.abs(th_now-this.th[i])<Math.PI)
			this.th[i] = (th_now + this.th[i]) * 0.5;
		else
			this.th[i] = (th_now + this.th[i]) * 0.5 + Math.PI;
		*/
		this.th[i] = th_now;
			
		
		// 回転行列の作成
		this.Re[i] = [[Math.cos(this.th[i]),-Math.sin(this.th[i])],[Math.sin(this.th[i]),Math.cos(this.th[i])]];
		

        // StiffnessWarpingの効果を無効にしたいときは
        // コメントアウトを外して回転行列を単位行列にする
        //this.Re[i] = numeric.identity(2);


		// 要素剛性マトリクス作成
		var KeTmp = numeric.rep([6,6],0);
		for(var bi=0; bi<3; bi++)
			for(var bj=0; bj<3; bj++)
				for(var j = 0; j < 2; j++) 
					for(var k = 0; k < 2; k++) 
						for(var l = 0; l < 2; l++) 
							for(var m = 0; m < 2; m++) 
								KeTmp[2 * bi + j][2 * bj + k] += this.Re[i][j][l]*this.Ke[i][2*bi+l][2*bj+m]*this.Re[i][k][m];


		// 全体剛性マトリクス作成
		for(var j=0; j<3; j++)
			for(var k=0; k<3; k++)
				for(var l=0; l<2; l++)
					for(var m=0; m<2; m++)
						this.K[2*this.tri[i][j]+l][2*this.tri[i][k]+m] += KeTmp[2*j+l][2*k+m];
		
		// 力オフセット作成のための回転剛性マトリクス作成
		var RKeTmp = numeric.rep([6,6],0);
		for(var bi = 0; bi < 3; bi++) 
			for(var bj = 0; bj < 3; bj++) 
				for(var j = 0; j < 2; j++) 
					for(var k = 0; k < 2; k++) 
						for(var l = 0; l < 2; l++) 
							RKeTmp[2*bi+j][2*bj+k] += this.Re[i][j][l] * this.Ke[i][2*bi+l][2*bj+k];

		for(var j=0; j<3; j++)
			for(var k=0; k<3; k++)
				for(var l=0; l<2; l++)
					for(var m=0; m<2; m++)
						RK[2*this.tri[i][j]+l][2*this.tri[i][k]+m] += RKeTmp[2*j+l][2*k+m];
						
	}

	
	// 力オフセット合成
	var initposVec = numeric.linspace(0,0,2*this.posNum);
	for(var i=0; i<this.posNum; i++)
		for(var j=0; j<2; j++)
			initposVec[2*i+j] = this.initpos[i][j];
	this.foffset = numeric.dot(RK,initposVec);
	this.foffset=numeric.neg(this.foffset);

}


// 境界条件を設定して変形計算を行う
// 境界条件は y=0 を固定，ノード番号spNodeに強制変位disp[2]を与える
// modeはStiffnessWarpingかNormal
FEM.prototype.calcDeformation = function(){
	
	// Stiffness Warping法の場合、剛性マトリクスを修正する
	this.makeMatrixKSW();
	
	var f = this.flist.length;
	var d = this.dlist.length;
	
	var Kff = numeric.rep([2*f,2*f],0);
	for(var i=0; i<f; i++)
		for(var j=0; j<f; j++)
			for(var k=0; k<2; k++)
				for(var l=0; l<2; l++)
					Kff[2*i+k][2*j+l] = this.K[2*this.flist[i]+k][2*this.flist[j]+l];
	
	var Kfd = numeric.rep([2*f,2*d],0);
	for(var i=0; i<f; i++)
		for(var j=0; j<d; j++)
			for(var k=0; k<2; k++)
				for(var l=0; l<2; l++)
					Kfd[2*i+k][2*j+l] = this.K[2*this.flist[i]+k][2*this.dlist[j]+l];

	var xd = numeric.linspace(0,0,2*d);
	for(var i=0; i<d; i++){
		for(var j=0; j<2; j++){
			xd[2*i+j] = this.initpos[this.dlist[i]][j] + this.ud[2*i+j];
		}
	}

	var felastic = numeric.linspace(0,0,2*f);
	for(var i=0; i<f; i++){
		for(var j=0; j<2; j++){
			felastic[2*i+j] = - this.foffset[2*this.flist[i]+j];
		}
	}
	
	var y = numeric.dot(Kfd,xd);
	y = numeric.sub(felastic,y);
	y = numeric.add(y, this.ff);
	var xf = numeric.solve(Kff,y);
	
	for(var i=0; i<f; i++)
		for(var j=0; j<2; j++)
			this.pos[this.flist[i]][j] = xf[2*i+j];
			
	for(var i=0; i<d; i++)
		for(var j=0; j<2; j++)
			this.pos[this.dlist[i]][j] = this.initpos[this.dlist[i]][j] + this.ud[2*i+j];
			
}


FEM.prototype.divideMatrices = function (dt) {

	var f = this.flist.length;
	var d = this.dlist.length;
			
	this.Kff = numeric.rep([2*f,2*f],0);
	for(var i=0; i<f; i++)
		for(var j=0; j<f; j++)
			for(var k=0; k<2; k++)
				for(var l=0; l<2; l++)
					this.Kff[2*i+k][2*j+l] = this.K[2*this.flist[i]+k][2*this.flist[j]+l];
	
	this.Kfd = numeric.rep([2*f,2*d],0);
	for(var i=0; i<f; i++)
		for(var j=0; j<d; j++)
			for(var k=0; k<2; k++)
				for(var l=0; l<2; l++)
					this.Kfd[2*i+k][2*j+l] = this.K[2*this.flist[i]+k][2*this.dlist[j]+l];

	this.Kdf = numeric.rep([2*d,2*f],0);
	for(var i=0; i<d; i++)
		for(var j=0; j<f; j++)
			for(var k=0; k<2; k++)
				for(var l=0; l<2; l++)
					this.Kdf[2*i+k][2*j+l] = this.K[2*this.dlist[i]+k][2*this.flist[j]+l];
	
	this.Kdd = numeric.rep([2*d,2*d],0);
	for(var i=0; i<d; i++)
		for(var j=0; j<d; j++)
			for(var k=0; k<2; k++)
				for(var l=0; l<2; l++)
					this.Kdd[2*i+k][2*j+l] = this.K[2*this.dlist[i]+k][2*this.dlist[j]+l];
	
	this.Mf = numeric.identity(2*f);
	for(var i=0; i<f; i++){
		this.Mf[2*i][2*i] = this.mass[this.flist[i]];
		this.Mf[2*i+1][2*i+1] = this.mass[this.flist[i]];
	}
		
	this.Md = numeric.identity(2*d);
	for(var i=0; i<d; i++){
		this.Md[2*i][2*i] = this.mass[this.dlist[i]];
		this.Md[2*i+1][2*i+1] = this.mass[this.dlist[i]];
	}
		
	// uf^{i}		
	this.uf = numeric.linspace(0,0,2*f);
	for(var i=0; i<f; i++)
		for(var j=0; j<2; j++)
			this.uf[2*i+j] = this.pos[this.flist[i]][j] - this.initpos[this.flist[i]][j];

	// vf^{i}		
	this.vf = numeric.linspace(0,0,2*f);
	for(var i=0; i<f; i++)
		for(var j=0; j<2; j++)
			this.vf[2*i+j] = this.Vel[2*this.flist[i]+j];

	// vd^{i}
	this.vdbuf = numeric.linspace(0,0,2*d);
	for(var i=0; i<d; i++)
		for(var j=0; j<2; j++)
			this.vdbuf[2*i+j] = this.Vel[2*this.dlist[i]+j];

	
	// vd^{i+1} (this.udはud^{i+1})
	this.vd = numeric.linspace(0,0,2*d);
	for(var i=0; i<d; i++)
		for(var j=0; j<2; j++)
			this.vd[2*i+j] = ( this.ud[2*i+j] + this.initpos[this.dlist[i]][j] - this.pos[this.dlist[i]][j] ) / dt;

	

	// xd^{i}
	this.xd = numeric.linspace(0,0,2*d);
	for(var i=0; i<d; i++)
		for(var j=0; j<2; j++)
			this.xd[2*i+j] = this.pos[this.dlist[i]][j];

	// xf^{i}	
	this.xf = numeric.linspace(0,0,2*f);
	for(var i=0; i<f; i++)
		for(var j=0; j<2; j++)
			this.xf[2*i+j] = this.pos[this.flist[i]][j];

	this.fof = numeric.linspace(0,0,2*f);
	for(var i=0; i<f; i++)
		for(var j=0; j<2; j++)
			this.fof[2*i+j] = this.foffset[2*this.flist[i]+j];

	this.fod = numeric.linspace(0,0,2*d);
	for(var i=0; i<d; i++)
		for(var j=0; j<2; j++)
			this.fod[2*i+j] = this.foffset[2*this.dlist[i]+j];

}


// 境界条件を設定して変形計算を行う
// 境界条件は y=0 を固定，ノード番号spNodeに強制変位disp[2]を与える
FEM.prototype.calcDynamicDeformation = function(dt){
	
	var f = this.flist.length;
	var d = this.dlist.length;

	this.makeMatrixKSW();
	this.divideMatrices(dt);
	
	var Mleft1 = numeric.mul(this.Mf,(1+this.alpha));
	var Mleft2 = numeric.mul(this.Kff,dt*(this.beta+dt));
	var Mleft = numeric.add(Mleft1,Mleft2);

	var Mright1 = numeric.dot(this.Mf,this.vf);
	var Mright2 = numeric.dot(this.Kff,this.xf);
	Mright2 = numeric.neg(Mright2);
	var tmp = numeric.dot(this.Kfd,this.xd);
	Mright2 = numeric.sub(Mright2,tmp);
	Mright2 = numeric.sub(Mright2,this.fof);
	Mright2 = numeric.add(Mright2,this.ff);
	Mright2 = numeric.mul(Mright2,dt);
	var Mright3 = numeric.mul(this.beta*dt+dt*dt, this.Kfd);
	tmp = numeric.dot(Mright3, this.vd);
	Mright3 = tmp;
	var Mright = numeric.add(Mright1,Mright2);
	Mright = numeric.add(Mright, Mright3);
	
	this.vf = numeric.solve(Mleft,Mright);

	for(var i=0; i<f; i++)
		for(var j=0; j<2; j++)
			this.Vel[2*this.flist[i]+j]=this.vf[2*i+j];

	for(var i=0; i<d; i++)
		for(var j=0; j<2; j++)
			this.Vel[2*this.dlist[i]+j]=(this.ud[2*i+j]-(this.pos[this.dlist[i]][j]-this.initpos[this.dlist[i]][j]))/dt;

	var duf = numeric.mul(dt,this.vf);
	this.uf = numeric.add(this.uf,duf);
	for(var i=0; i<f; i++)
		for(var j=0; j<2; j++)
			this.pos[this.flist[i]][j] = this.initpos[this.flist[i]][j] + this.uf[2*i+j];
	
	for(var i=0; i<d; i++)
		for(var j=0; j<2; j++)
			this.pos[this.dlist[i]][j] = this.initpos[this.dlist[i]][j] + this.ud[2*i+j];

	// 外力の計算
	var tmp1, tmp2, tmp3, tmp4, tmp5;
	var tmp23, tmp1234;
	var tmpDot;
	// 速度に関する項を外力計算に入れると値が
	// 大きすぎる値になったのでコメントアウトした
	/*
	tmp1 = numeric.mul(dt*this.beta+dt*dt, this.Kdf);
	tmpDot = numeric.dot(tmp, this.vf);
	tmp1 = numeric.clone(tmpDot);

	tmp2 = numeric.mul(1+dt*this.alpha, this.Md);

	tmp3 = numeric.mul(dt*this.beta+dt*dt, this.Kdd);

	tmp23 = numeric.add(tmp2, tmp3);
	tmpDot = numeric.dot(tmp23, this.vd);
	tmp23 = numeric.clone(tmpDot);

	tmp4 = numeric.dot(this.Md, this.vdbuf);
	tmp4 = numeric.neg(tmp4);

	tmp1234 = numeric.add(tmp1, tmp23);
	tmp1234 = numeric.add(tmp1234, tmp4);
	tmp1234 = numeric.mul(tmp1234, 1/dt);
	*/

	tmp5 = numeric.dot(this.Kdf, this.xf);
	tmp6 = numeric.dot(this.Kdd, this.xd);

	var fd = numeric.add(tmp5, tmp6);

	// 速度を考慮する場合は以下のコメントアウトを外す
	// fd = numeric.add(fd, tmp1234);

	fd = numeric.add(fd, this.fod);


	for(var i=0; i<d; i++){
		this.f[2*this.dlist[i]+0] = fd[2*i+0];
		this.f[2*this.dlist[i]+1] = fd[2*i+1];
	}
}	


FEM.prototype.modifyPosCld = function(xmin, ymin, xmax, ymax){
	for(var i=0; i<this.pos.length; i++) {
		if(this.pos[i][0]<xmin) {
			this.pos[i][0]=xmin;
			this.Vel[2*i] = 0;
			this.Vel[2*i+1] = 0;
		}
		if(this.pos[i][0]>xmax) {
			this.pos[i][0]=xmax;
			this.Vel[2*i] = 0;
			this.Vel[2*i+1] = 0;
		}
		if(this.pos[i][1]<ymin) {
			this.pos[i][1]=ymin;
			this.Vel[2*i] = 0;
			this.Vel[2*i+1] = 0;
		}
		if(this.pos[i][1]>ymax) {
			this.pos[i][1]=ymax;
			this.Vel[2*i] = 0;
			this.Vel[2*i+1] = 0;
		}
	}
}


FEM.prototype.getForce = function () {
	var force = Array(this.holdNode.length);
	for(var i=0; i<this.holdNode.length; i++){
		force[i] = [0,0];
		for(var j = 0; j < this.holdNode[i].length; j++) {
			force[i][0] += this.f[2*this.holdNode[i][j]+0];
			force[i][1] += this.f[2*this.holdNode[i][j]+1];
		}
	}
	return force;
}
	

FEM.prototype.calcStress = function () {
    for (var i = 0; i < this.tri.length; i++) {
        var xe = [0,0,0,0,0,0];
        for(var j=0; j<3; j++){
            xe[2*j] = this.pos[this.tri[i][j]][0];
            xe[2*j+1] = this.pos[this.tri[i][j]][1];
        }
        var x0e = [0,0,0,0,0,0];
        for(var j=0; j<3; j++){
            x0e[2*j] = this.initpos[this.tri[i][j]][0];
            x0e[2*j+1] = this.initpos[this.tri[i][j]][1];
        }
        //var ReInv = numeric.transpose(this.Re[i]);
		var ReInv = numeric.rep([6,6],0);
		for(var j=0; j<3; j++)
			for(var k=0; k<2; k++)
				for(var l = 0; l < 2; l++) 
					ReInv[2*j+k][2*j+l] = this.Re[i][l][k];
        var strain = numeric.dot(ReInv,xe);
        strain = numeric.sub(strain,x0e);
        var tmp = numeric.dot(this.De[i],this.Be[i]);
        var stress = numeric.dot(tmp,strain);
        var sigma1 = (stress[0]+stress[1])*0.5+Math.sqrt((stress[0]-stress[1])*(stress[0]-stress[1])*0.25+stress[2]*stress[2]);
        var sigma2 = (stress[0]+stress[1])*0.5-Math.sqrt((stress[0]-stress[1])*(stress[0]-stress[1])*0.25+stress[2]*stress[2]);
        if(Math.abs(sigma1)>Math.abs(sigma2)) {
            this.maxPStress[i]=Math.abs(sigma1);
        } else {
            this.maxPStress[i]=Math.abs(sigma2);
        }
    }
}


FEM.prototype.removeElement=function () {
    for(var i=0; i<this.tri.length; i++) {
        if(this.maxPStress[i]>this.thrPStress) {
            this.removedFlag[i] = true;
        }
    }
}


// 三角形定ひずみ要素の要素Bマトリクスを作成する関数
function makeMatrixB(p1,p2,p3){
	var Be = [[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]];
	var mat = [[1,p1[0],p1[1]], [1,p2[0],p2[1]], [1,p3[0],p3[1]]];
	var delta = numeric.det(mat);
	var dd = 1.0/delta;
	Be[0][0] = (p2[1]-p3[1])*dd;
	Be[0][2] = (p3[1]-p1[1])*dd;
	Be[0][4] = (p1[1]-p2[1])*dd;
	Be[1][1] = (p3[0]-p2[0])*dd;
	Be[1][3] = (p1[0]-p3[0])*dd;
	Be[1][5] = (p2[0]-p1[0])*dd;
	Be[2][0] = (p3[0]-p2[0])*dd;
	Be[2][1] = (p2[1]-p3[1])*dd;
	Be[2][2] = (p1[0]-p3[0])*dd;
	Be[2][3] = (p3[1]-p1[1])*dd;
	Be[2][4] = (p2[0]-p1[0])*dd;
	Be[2][5] = (p1[1]-p2[1])*dd;
	return Be;
}
