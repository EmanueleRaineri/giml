type segment = { pos:int array; score:float array; mean_score:float; dist:float};;


let sum v =  
	Array.fold_left (+.) 0. v
;;

let mean v =
	let l = float_of_int (Array.length v)
	in (sum v) /. l
;;

let map2 f v1 v2 =
	let rec aux_map2 v1 v2 acc idx =
		if (idx=Array.length v1)
		then acc
		else 
		aux_map2 v1 v2 (acc +. (f v1.(idx) v2.(idx))) idx
	in
	aux_map2 v1 v2 0. 0
;;

let diff v1 v2 =
	let l1 = Array.length v1 and 
	l2 = Array.length v2 in
	ignore l2;
	let res = Array.init l1 (fun i -> 0.) in
	for i = 0 to l1 -1 do
		res.(i) <- ( v1.(i) -. v2.(i) )
	done;
	res
;;

let prod v1 v2 =
	let l1 = Array.length v1 and 
	l2 = Array.length v2 in
	ignore l2;
	let res = Array.init l1 (fun i -> 0.) in
	for i = 0 to l1 -1 do
		res.(i) <- ( v1.(i) *. v2.(i) )
	done;
	res
;;

let norm v =
	let p = prod v v in
	sum p
;;


let make_segment pos score =
	let mean_score = mean score	
	in let v = Array.map (fun x -> x -. mean_score) score
	in 
	let dist = norm v in
	{ pos=pos; score=score; mean_score=mean_score; dist = dist  }
;;
