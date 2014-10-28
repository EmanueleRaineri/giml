exception Exit;;

let parent i = i/2;;

let left i = 2*i;;

let right i = 2*i+1;;

let   max_heapify a i hsize =
	let idx = ref i in
	let largest = ref (-1) in
	try
	while ( true ) do
		let l = left !idx 
		and r = right !idx 
		in
		if (l<=hsize && a.(l-1)>a.(!idx -1 )) then
			largest:=l
		else
			largest:=!idx;
		if (r<=hsize && a.(r-1)>a.(!largest-1)) then
			largest:=r;
		if (!largest != !idx) then
			let tmp = a.(!idx-1)
			in begin 
				a.(!idx-1) <- a.(!largest-1);
				a.(!largest-1)<-tmp;
				idx := !largest
			end
		else
			raise Exit
	done;
	a
	with Exit -> a
;;

let b = [|16;4;10;14;7;9;3;2;8;1|] and hsize=10 in
	let h = max_heapify b 2 hsize 
	in Array.iter (Printf.printf " %d ") h;;
