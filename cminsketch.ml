(* Multiply shift is done on normal ints.  These are not necessarily
 * 32-bit or 64-bit (usually, they have one bit less precision), so
 * care must be taken. *)

(* Warning: If this value is too small you *will* get out of bounds
 * errors.  Thus, a mythical OCaml implementation that does not
 * reserve any bits for bookkeeping will not handle this properly.
 * Fortunately, this is unlikely to change. *)
let int_size = Sys.word_size - 1
let multiply_shift m a x = (a * x) lsr (int_size - m)

let euler = exp 1.
(* XXX probably there is something more efficient, but this is only
 * called for sketch creation, which is fairly infrequent *)
let lg x = log x /. log 2.
let int_ceil x = int_of_float (ceil x)
(* Sys.word_size - 2 bits is "just right" (30 bits for 32-bit
 * and 62 bits for 64-bit), because one bit is thrown out due to
 * implementation details, and one bit is thrown out due to us
 * specifying the number to be odd *)
let random_odd_int () =
    if Sys.word_size = 32
    then Random.bits () * 2 + 1 (* nice and simple *)
    (* This implementation "wastes" 28 bits, but IIUC there's
     * no way to ask for less entropy, since all of the functions
     * in the Random module use bits internally *)
    else (Random.bits () lor (Random.bits () lsl 30) lor (Random.bits () lsl 60)) * 2 + 1
let step_matrix a i j c = a.(i).(j) <- a.(i).(j) + c
let minimum a = Array.fold_left min max_int a

type cminsketch = CMinSketch of int * int array array * int array

(* utilizes RNG *)
(* arguments are to be interpeted as thus: the error in answering
 * a query is within a factor of epsilon with probability delta.
 * You get more accurate results for small epsilon and large delta.
 * Of course, don't set epsilon = 0 or delta = 1: then you're
 * no longer using a probabilistic data structure. *)
let make epsilon delta =
    (* We fudge the width to be a little larger, ensuring that it
     * is a power of two for the benefit of our algorithm.  This means
     * the actual epsilon you will get is smaller than what you
     * originally specified. *)
    (* XXX Add error checking for the parameters *)
    let m = int_ceil (lg (euler /. epsilon)) in
    let width = 1 lsl m
    and depth = int_ceil (log (1. /. delta)) in
    CMinSketch
        (m,
         Array.make_matrix depth width 0,
         Array.init depth (fun _ -> random_odd_int ()))

let update (CMinSketch (m, sketch, hfs)) ix c =
    Array.iteri (fun i a -> step_matrix sketch i (multiply_shift m a ix) c) hfs

let query (CMinSketch (m, sketch, hfs)) ix =
    (* No fusion :-( so this generates an intermediate data structure
     * that is immediately discarded.  We could probably write a
     * minimum_mapi function but meh *)
    minimum (Array.mapi (fun i a -> sketch.(i).(multiply_shift m a ix)) hfs)


(* To implement: *)
(* range_query, inner_product_query *)
(* phi_quantiles, heavy_hitters *)

let () =
    let x = make 0.998 0.002 in
    update x 3 4;
    update x 3 2;
    update x 24435 1;
    update x 2323434 1;
    update x 223434 1;
    print_int (query x 223434);
    print_string "\n";
    ()
