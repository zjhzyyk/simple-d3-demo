/**
math 2605 project: Jacobi Algorithm.
author: Yukang Yang
*/

// alert("hello");
var max = 500;
var min = -500;
var Matrix = function(){
	this.width=0;
	this.height=0;
	this.entries=[];
};

function point(a,b){
    this.x = a;
    this.y = b;
}

Matrix.prototype = {
	/**
	initialize matrix instance with given dimension
	**/
	setDimension : function(w,h) {
		this.entries = [];
		for (var i=0; i<h; i++)
			this.entries.push(new Array(w));
		this.width = w;
		this.height = h;
	},
	/**
	clone to target matrix
	*/
	copyTo : function(target) {
		target.setDimension(this.width, this.height);
		for (var i = 0; i<this.height; i++) for (var j = 0; j<this.width; j++) target.entries[i][j] = this.entries[i][j];
	},
	/**
	randomly generate a 5 by 5 matrix
	*/
	generate : function(){
		this.setDimension(5,5);
		for (var i=0; i<5; i++) {
			for (var j=0; j<5; j++) {
				if (j<i) this.entries[i][j] = this.entries[j][i]; //for the sake of symmetric matrix
				else this.entries[i][j] = Math.floor(Math.random()*(max-min+1))+min; //generate an integer in [min, max]
			}
		}
	},
	/**
	minus current matrix by m
	*/
	minus : function (m) {
		var ret = new Matrix();
		this.copyTo(ret);
		if (this.width==m.width && this.height==m.height) {
			for (var i = 0; i<this.height; i++) for (var j = 0; j<this.width; j++) ret.entries[i][j] = this.entries[i][j]-m.entries[i][j];
		}
		return ret;
	},
	scalarMultiply : function (scalar) {
		var ret = new Matrix();
		this.copyTo(ret);
		for (var i = 0; i<this.height; i++) for (var j = 0; j<this.width; j++) ret.entries[i][j] *= scalar;
		return ret;
	},
	/**
	precondition: i<j
	orthgonally diagonalize a 2 by 2 submatrix in a symmetric matrix and U.
	*/
	diag : function(i, j){
		var a = this.entries[i][i];
		var b = this.entries[i][j];
		var d = this.entries[j][j];
		var mat = new Matrix();
		mat.setDimension(2,2);
		mat.entries[0][0] = a;
		mat.entries[1][0] = mat.entries[0][1] = b;
		mat.entries[1][1] = d;
		var idt = new Matrix();
		idt.setDimension(2,2);
		idt.entries[0][0] = idt.entries[1][1] = 1;
		idt.entries[1][0] = idt.entries[0][1] = 0;
		//calc the larger eigenvalue
		var ev = (a+d)/2.0 + Math.sqrt(b*b+Math.pow((a-d)/2.0,2));
		mat = mat.minus(idt.scalarMultiply(ev));
		mat.entries[0][1] = mat.entries[0][0];
		mat.entries[1][1] = mat.entries[1][0];
		mat.entries[0][0] = mat.entries[1][1];
		mat.entries[1][0] = -mat.entries[0][1];
		mat = mat.scalarMultiply(1.0/Math.sqrt(mat.entries[0][0]*mat.entries[0][0]+mat.entries[0][1]*mat.entries[0][1]));
		return mat;
	},
	transpose : function(){
		var ret = new Matrix();
		ret.setDimension(this.height, this.width);
		for (var i = 0; i<this.height; i++) for (var j = 0; j<this.width; j++) ret.entries[j][i] = this.entries[i][j];
		return ret;
	},
	matrixMultiply : function(b){
		if (b.height != this.width) return;
		var ret = new Matrix();
		ret.setDimension(this.height, b.width);
		for (var i = 0; i < this.height; i++) for (var j = 0; j<b.width; j++) {
			ret.entries[i][j] = 0;
			for (var k = 0; k<this.width; k++) {
				ret.entries[i][j]+= (this.entries[i][k]*b.entries[k][j]);
			}
		}
		return ret;
	},
	/**
	calc off value
	*/
	off : function() {
		var ret = 0;
		for (var i = 0; i<this.height; i++) for (var j =0; j<this.width; j++)
			if (i!=j) ret+= (this.entries[i][j]*this.entries[i][j]);
		return ret;
	},
	getRotMat : function(u,a,b) {
		var ret = new Matrix();
		ret.setDimension(this.height, this.height);
		for (var i = 0; i<this.height; i++) for (var j = 0; j<this.height; j++) {
			if (i==j) ret.entries[i][j] = 1;
			else ret.entries[i][j] = 0;
		}
		ret.entries[a][a] = u.entries[0][0];
		ret.entries[a][b] = u.entries[0][1];
		ret.entries[b][a] = u.entries[1][0];
		ret.entries[b][b] = u.entries[1][1];
		return ret;
	},
	jacobi : function(sorting){
		var ret = new Matrix();
		var offArr = [];
		var tmp=100.0;
		var i = 0, count = 3000;
		this.copyTo(ret);
		var mi=0;
		var mj=ret.width;
		while (tmp>1e-9 && count--) {
			var max = -10000;
			if (sorting === true) {
			//find largest off-diagonal entry
				for (i= 0; i<ret.height; i++) for (var j = i+1; j<ret.width; j++)
					if (Math.abs(ret.entries[i][j])>max) {
						mi = i;
						mj = j;
						max = Math.abs(ret.entries[i][j]);
					}
			}
			else {
				if (mj-1==mi) {
					mj = ret.width-1;
					if (mi+1==4)
						mi = 0;
					else
						mi++;
				}
				else
					mj--;
				console.log("mi: " + mi.toString() + " mj: " + mj.toString());
			}
			var g = ret.getRotMat(ret.diag(mi,mj), mi, mj);
			ret = g.transpose().matrixMultiply(ret).matrixMultiply(g);
			// console.log(JSON.stringify(ret.entries));
			tmp = ret.off();
			offArr.push(tmp);
		}
		if (count===-1 || count===0) {
			$("eigenvalues").innerHTML = "fail";
			return;
		}
		var val = "";
		for (i=0; i<ret.width; i++) val+=ret.entries[i][i].toString()+"<br />";
		$("eigenvalues").innerHTML = val;
		return offArr;
	}
};

function $(id) {
	return document.getElementById(id);
}

var gwidth = 750,
	gheight = 500,
	xmin = -20,
	ymin = -20;

function init(){
	var omtTable = $("omat");
	while (omtTable.rows.length!=5) {
		row = omtTable.insertRow(0);
		for (i=0; i<5; i++) row.insertCell(0);
	}
	//create graph for current matrix
	var chart = d3.select("body")
		.append("svg:svg")
		.attr("width", gwidth)
		.attr("height",gheight)
		.attr("id", "cur_graph");
	//create g
	var g = chart.append("svg:g")
		.attr("transform", "translate(" + (-xmin).toString() + "," + (gheight+ymin).toString()+") matrix(1 0 0 -1 0 0)");
}

/**
main part
*/
function main(){
	var sorting = document.querySelector("#optionsRadios1").checked;
	console.log("sorting: " + sorting.toString());
	var mat = new Matrix();
	mat.generate();
	var omt = "";
	var i=0;
	var j=0;
	var omtTable = $("omat");
	var row;
	var cell;
	for (i=0;i<5; i++) {
		cell = omtTable.rows[i].cells;
		for (j=0;j<5;j++) {
			cell[j].align = "center";
			cell[j].valign = "middle";
			cell[j].innerHTML = mat.entries[i][j];
		}
	}
	// mat.entries = [
	// 	[309,-114,-896,-16,255],
	// 	[-114,8,199,-62,-695],
	// 	[-896,199,-449,-153,579],
	// 	[-16,-62,-153,-961,382],
	// 	[255,-695,579,382,294]
	// ];
	// mat.width = 5;
	// mat.height = 5;
	var off = mat.jacobi(sorting);
	//draw the graph
	var k = Math.log(9/10),
		b = Math.log(mat.off());
	// console.log("k"+k.toString());
	// console.log("b"+b.toString());
	var fx = function(x) {
			return k*x+b;
		},
		fxinverse = function(y) {
			return (y-b)/k;
		};
	var scale = 21/20;
	var xl = fxinverse(fx(0)*scale),
		xr = fxinverse(0)*scale,
		yl = fx(xr),
		yr = fx(0)*scale;
	// console.log("xr"+xr.toString());
	var xscale = d3.scale.linear().domain([xl, xr]).range([xl*(gwidth+xmin)/xr, gwidth+xmin]),
		yscale = d3.scale.linear().domain([yl, yr]).range([yl*(gheight+ymin)/yr, gheight+ymin]);
	var g = document.querySelector("#cur_graph g");
	//first clear the graph of last matrix
	while (g.firstChild!==null) {
		g.removeChild(g.firstChild);
	}
	g = d3.select("#cur_graph g");
	//create axis
	//x
	g.append("svg:line")
		.attr("x1", xmin)
		.attr("y1", 0)
		.attr("x2", gwidth+xmin)
		.attr("y2", 0)
		.attr("class", "axis");
	//y
	g.append("svg:line")
		.attr("x1", 0)
		.attr("y1", ymin)
		.attr("x2", 0)
		.attr("y2", gheight+ymin)
		.attr("class", "axis");
	//draw 12 axis labels
	function xtick(){
		var arr = [];
		var interval = Math.floor(fxinverse(0)/12);
		var i = interval;
		for (;i<fxinverse(0); i+=interval) {
			arr.push(i);
		}
		return arr;
	}
	function ytick(){
		var arr = [];
		var interval = Math.floor(fx(0)/10);
		var i = interval;
		for (;i<fx(0); i+=interval) {
			arr.push(i);
		}
		return arr;
	}
	g.selectAll(".xLabel")
		.data(xtick())
		.enter().append("svg:text")
		.attr("transform", "matrix(1 0 0 -1 0 0)")
		.attr("class", "xLabel")
		.text(function(d){return d.toString();})
		.attr("x", function(d) { return xscale(d); })
		.attr("y", 15)
		.attr("text-anchor", "middle");

	g.selectAll(".yLabel")
		.data(ytick())
		.enter().append("svg:text")
		.attr("transform", "matrix(1 0 0 -1 0 0)")
		.attr("class", "yLabel")
		.text(function(d){return d.toString();})
		.attr("x", 0)
		.attr("y", function(d) { return -yscale(d); })
		.attr("text-anchor", "right")
		.attr("dx", -18)
		.attr("dy", 4);

	g.selectAll(".xTicks")
		.data(xtick())
		.enter().append("svg:line")
		.attr("class", "xTicks")
		.attr("x1", function(d) { return xscale(d); })
		.attr("y1", 0)
		.attr("x2", function(d) { return xscale(d); })
		.attr("y2", 6);

	g.selectAll(".yTicks")
		.data(ytick())
		.enter().append("svg:line")
		.attr("class", "yTicks")
		.attr("y1", function(d) { return yscale(d); })
		.attr("x1", 6)
		.attr("y2", function(d) { return yscale(d); })
		.attr("x2", 0);
	//draw function, i.e. the upper bound of off value
	g.append("svg:line")
		.attr("x1", xscale(xl))
		.attr("y1", yscale(yr))
		.attr("x2", xscale(xr))
		.attr("y2", yscale(yl))
		.attr("class", "function");
	//draw off value
	var l = off.length;
	console.log(off);
	for (i = 0; i<l; i++) {
		g.append("svg:circle")
			.attr("cx", xscale(i+1))
			.attr("cy", yscale(Math.log(off[i])))
			.attr("r", 5);
	}
}
