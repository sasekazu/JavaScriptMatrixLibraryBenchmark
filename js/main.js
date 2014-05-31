// JavaScript Document
/// <reference path="jquery.js" />
/// <reference path="numeric-1.2.6.min.js" />


$(document).ready(function () {

	var init = 5;
	var max = 500;
	var bin = 20;
	var data = [];
	for(var N = init; N <= max ; N += bin) {
		var A = numeric.random([N, N]);
		var b = numeric.random([N]);
		var time0 = new Date();
		numeric.solve(A, b);
		var time1 = new Date();
		console.log("N=" + N + " " + (time1 - time0) + " [ms]");
		data.push([N, time1-time0]);

		/*
		// 計算結果を確認したいときは以下の3行を有効にする
		var x = numeric.solve(A, b);
		var AxMinusb = numeric.sub(numeric.dot(A,x),b);
		console.log("A*x-b = " + AxMinusb);
		*/
	}

	// グラフ作成
	$("#canvas-frame").empty();
	var plotDatas = [];
	plotDatas[0] = data;
	var options = {
		axes: {
			xaxis: {
				label: 'N',
				min: 0,
				max: data[data.length - 1][0],
				tickInterval: bin * 10
			},
			yaxis: {
				label: 'Time [ms]',
				min: 0
			}
		},
		seriesDefaults: {
			showMarker: true
		}
	};
	$.jqplot("canvas-frame", plotDatas, options);

	
} );

