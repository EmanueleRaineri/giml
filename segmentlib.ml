type segment = { pos:int array; score:float array; mean_score:float; dist:float};;

let print_float_array v =
 let n = Array.length v in
 for i = 0 to (n-2) do
 	Printf.fprintf stdout "%d:%.4f\t" i v.(i)
 done;
 	Printf.fprintf stdout "%d:%.4f\n" (n-1) v.(n-1)
;;

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
		aux_map2 v1 v2 (acc +. (f v1.(idx) v2.(idx))) (idx + 1)
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

let l2dist score =
	let mean_score = mean score in
	let v = Array.map (fun x -> x -. mean_score) score
	in 
	 norm v
;;

let make_segment pos score =
	let mean_score = mean score	
	in let v = Array.map (fun x -> x -. mean_score) score
	in 
	let dist = norm v in
	{ pos=pos; score=score; mean_score=mean_score; dist = dist  }
;;

let join_vectors v1 v2 =
	let l1 = Array.length v1 
	and l2 = Array.length v2
	in let jv =  Array.init ( l1 + l2 ) 
	(fun i -> if ( i < l1 ) then v1.(i) else v2.(i - l1 ))
	in
	jv
;;

let join_segments seg1 seg2 =
	 let jpos = join_vectors seg1.pos seg2.pos 
	and jscore = join_vectors seg1.score seg2.score in 
	let jms = mean jscore
	and jdist = l2dist jscore 
	in
	{ pos=jpos; score=jscore; mean_score=jms; dist=jdist   }
;;

let delta_segments seg1 seg2 lambda =
	let jsco = join_vectors seg1.score seg2.score in
	let joinedl2dist = l2dist jsco in
	joinedl2dist -. seg1.dist -. seg2.dist -. lambda	
;;

let range a b =
	let l = b-a+1
	in Array.init l (fun i -> (a + i ))
;;

let all_pairs segv lambda =
	let les = Array.length segv in
	let cmp = fun  i -> 
	let seg1 = segv.(i) 
	and seg2=segv.(i+1)
	in  delta_segments seg1 seg2 lambda
	in
	Array.map cmp (range 0 (les-2))
;;

let maxi v =
	let l = Array.length v in
	let idx = ref 0 in
	for i=0 to (l-1) do
		if v.(i)>=v.(!idx) then idx:=i
	done;
	!idx
;;

let mini v =
	let l = Array.length v in
	let idx = ref 0 in
	for i=0 to (l-1) do
		if v.(i)<=v.(!idx) then idx:=i
	done;
	!idx
;;

let update_segmentation segv ap =
	let idx = mini ap
	in
	if (ap.(idx)<0) 
	then begin
		let new_segv = Array.init (Array.length ap -1 )
		(fun i -> segv.(i))
		in 
		new_segv.(idx) <- join_segments segv.(idx) segv.(idx+1);
		for i =idx+1 to Array.length ap -2 do
			new_segv.(i)<-segv.(i)
		done;
		new_segv
	end
	else
		segv
;;

let print_segment seg =
	let min_idx = seg.pos.(mini(seg.pos))
	and max_idx = seg.pos.(maxi(seg.pos)) in
	Printf.fprintf stdout  "%d\t%d\t" min_idx max_idx;
	Printf.fprintf stdout  "%.4f\t" seg.mean_score;
	Printf.fprintf stdout "%.4f\n" seg.dist
;;
