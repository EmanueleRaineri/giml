open Segmentlib;;


let _ =
	Printf.printf "ciao!\n";
	let seg1=  Segmentlib.make_segment [|10|] [|34.|]
	and seg2 = Segmentlib.make_segment [|20|] [|34.|]
	and seg3 = Segmentlib.make_segment [|30|] [|15.|]
	in
	let seg_array = [|seg1;seg2;seg3|] and 
	lambda = 5.0 in
	let ap =  Segmentlib.all_pairs seg_array lambda
	in
	Segmentlib.print_float_array ap;
	Segmentlib.print_segment seg1;
	Array.iter (fun s -> Segmentlib.print_segment s) seg_array
;;
