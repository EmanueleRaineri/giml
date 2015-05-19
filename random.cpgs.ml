open Random;;

let split_tab = Str.split (Str.regexp "\t")
;;

let  upload inch   =
			let acc = ref [] 
			and pos = ref [] in	
			let rec aux_upload () =
			let line = 
				input_line inch
				in
        		let fields = Array.of_list ( split_tab line ) in 
				let coord = int_of_string(fields.(1))
        		and nc     = float_of_string(fields.(5))
        		and c      = float_of_string(fields.(6))
				in  
					acc := (nc /. (nc +. c))::!acc  ;
					pos := coord::!pos;
					aux_upload ()
		in 
		try
		aux_upload ()
		with End_of_file ->  (Array.of_list (List.rev !pos), 
			Array.of_list (List.rev !acc))
;;

let avg a =
	let ub = Array.length a
	and sum = ref 0. in
	 for i = 0 to ub -1 do
		sum := !sum +. a.(i)
	done;
	!sum /. (float_of_int ub)
;;

let var a m =
	let ub = Array.length a
	and sum = ref 0. in
	 for i = 0 to ub -1 do
		sum := !sum +. ( a.(i) -. m )**2.0
	done;
	!sum /. (float_of_int ub)
;;


let _ =
	let le =15 in
        let (pos, met) = upload stdin 
		in   let bound = Array.length met 
		in 
		(*Printf.printf "%d\n" bound;*)
		for i = 0 to 10000 do
			let irnd1 = Random.int (bound - 100) in
			let slice_met = Array.sub met irnd1 le
			and slice_pos = Array.sub pos irnd1 le
			in let  m = avg slice_met
			in let  v = var slice_met m
			in Printf.printf "%d\t%d\t%f\t%g\n" ( slice_pos.(0) ) ( slice_pos.(le-1) ) m v
		done
;;

