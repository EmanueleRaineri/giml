let split_comma = Str.split(Str.regexp ",");;

let split_tab   = Str.split (Str.regexp "\t");;

let parse_values s =
	let a = Array.of_list (split_comma s )
	in Array.map (fun e -> int_of_string e) a
;;

let gammaln ( z : float ) =
    let cof =[|76.18009172947146;-86.50532032941677;
              24.01409824083091;-1.231739572450155;
              0.1208650973866179e-2;-0.5395239384953e-5|]
    in let y = ref z in
    let x = z in
    let tmp=x +.5.5 in
    let tmp = tmp -. (x +. 0.5)*.log(tmp) in
    let ser= ref 1.000000000190015 in
    for j= 0 to 5  do 
        y:=!y +. 1.0;
        ser := !ser +. cof.(j) /. !y 
    done;
-. tmp +. log(2.5066282746310005*. !ser /. x) 
;;

let factln n  =
     gammaln (float_of_int n +. 1.0)
;;

let logbico n k   = 
     factln n -. factln k -. factln ( n - k) 
;;

let logpow ( e:float ) ( base:float ) = 
	if ( e = 0.0 && base = 0.0 ) then 0.0 
	else ( e *. ( log base ) )  
;;

(* binomial probability  (n choose k) p^k q^(n-k) *)

let logpbico (n:int) (k:int) (p:float)  =
    let q = 1.0 -. p in
	logbico n k +. 
    logpow ( float_of_int k ) p +. 
    logpow ( float_of_int n -. float_of_int k ) q
;;

let range l =
	let v = Array.make  l 0 in 
	for i = 0 to (l-1) do
		v.(i)<-i
	done;
	v
;;

let loglikarray nc c theta=
	let l = Array.length nc
	in let r = range l
	in Array.map
	(fun i -> logpbico (nc.(i) + c.(i)) nc.(i) theta) r 
;;

let _ =
	try 
	while true do
		let s = input_line stdin
		in let fields = Array.of_list (split_tab s) 
		in (*Printf.printf "%d\n" (Array.length fields)*)
		let nc1 = parse_values fields.(3)
		and c1  = parse_values fields.(4)
		and nc2 = parse_values fields.(5)
		and c2  = parse_values fields.(6)
		in
		let sumnc1 = Array.fold_left (+) 0 nc1
		and sumc1 = Array.fold_left (+) 0 c1
		and sumnc2 = Array.fold_left (+) 0 nc2
		and sumc2 = Array.fold_left (+) 0 c2
		in
		let fsumnc1 = float_of_int sumnc1
		and fsumc1  = float_of_int sumc1
		and fsumnc2 = float_of_int sumnc2
		and fsumc2  = float_of_int sumc2
		in let theta1 = fsumnc1 /. ( fsumnc1 +. fsumc1 )
		and theta2 = fsumnc2 /. ( fsumnc2 +. fsumc2 )
		and theta12 = ( fsumnc1 +. fsumnc2 ) /. 
			( fsumnc1 +. fsumc1 +. fsumnc2 +.  fsumc2 ) 
		in
		let l1 = loglikarray nc1 c1 theta1
		and l2 = loglikarray nc2 c2 theta2
		in
		let suml1 = Array.fold_left ( +. ) 0. l1
		and suml2 = Array.fold_left ( +. ) 0. l2 in
		let nc12 = Array.append nc1 nc2 
		and c12 = Array.append c1 c2 
		in let l12 = loglikarray nc12 c12 theta12
		in let suml12 = Array.fold_left ( +. ) 0. l12 in
		Printf.printf "%s\t%s\t%s\t" fields.(0) fields.(1) fields.(2);
		Printf.printf "%d\t" ( Array.length nc1 );
		Printf.printf "%d\t" ( Array.length c1 );
		Printf.printf "%d\t" ( Array.length nc2 );
		Printf.printf "%d\t" ( Array.length c2 );
		Printf.printf "%f\t%f\t%f\t" theta1 theta2 theta12;
		Printf.printf "%f\t%f\t" suml1 suml2;
		Printf.printf "%f\t" suml12;
		Printf.printf "%f\n" ((suml1 +. suml2) -. suml12);
	done;
	()
	with End_of_file ->()
;;
