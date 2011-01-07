(* Multiply shift is done on normal ints.  These are not necessarily
 * 32-bit or 64-bit (usually, they have one bit less precision), so
 * care must be taken. *)

(** Number of bits in an OCaml int (not native int.) *)
(* Warning: If this value is too small you *will* get out of bounds
 * errors.  Thus, a mythical OCaml implementation that does not
 * reserve any bits for bookkeeping will not handle this properly.
 * Fortunately, this is unlikely to change. *)
let int_size = Sys.word_size - 1
(* The original paper suggested using Carter and Wegman's universal hash
 * family, which was strongly 2-independent, and the original paper
 * stipulates this requirement.  However, it has only proven that
 * multiply shift is weakly universal.  Fortunately, the proofs
 * in the paper don't require 2-independence, and the author of the
 * paper verified that multiply shift should have the necessary
 * theoretic properties. *)
let multiply_shift ~m ~a ~x = (a * x) lsr (int_size - m)

(** Euler's constant *)
let euler = exp 1.

(** Base 2 logarithm *)
(* XXX probably there is something more efficient, but this is only
 * called for sketch creation, which is fairly infrequent *)
let lg x = log x /. log 2.

(** Rounds a float up to the nearest integer. *)
let int_ceil x = int_of_float (ceil x)

(** Generates a random odd int.  May be negative. *)
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

(** Increases a the value at [(i, j)] in matrix [a] by [c]. *)
let step_matrix a i j c = a.(i).(j) <- a.(i).(j) + c

(** Finds the minimum value of an array. *)
let minimum a = Array.fold_left min max_int a

type sketch = { lg_width : int;
                count: int array array;
                hash_functions : int array }

let make ~epsilon ~delta =
    if epsilon <= 0.0 then invalid_arg "Cminsketch.make: epsilon must be greater than 0.0" else
    if delta >= 1.0 then invalid_arg "Cminsketch.make: delta must be less than 1.0" else
    if delta <= 0.0 then invalid_arg "Cminsketch.make: delta must be greater than 0.0" else
    (* We fudge the width to be a little larger, ensuring that it
     * is a power of two for the benefit of our algorithm.  This means
     * the actual epsilon you will get is smaller than what you
     * originally specified. *)
    let m = int_ceil (lg (euler /. epsilon)) in
    if m < 0 then failwith "Cminsketch.make: internal error, lg_width less than 0" else
    let width = 1 lsl m
    and depth = int_ceil (log (1. /. delta)) in
    if width <= 0 then failwith "Cminsketch.make: internal error, width less than 1" else
    if depth <= 0 then failwith "Cminsketch.make: internal error, depth less than 1" else
    { lg_width = m;
      count = Array.make_matrix depth width 0;
      hash_functions = Array.init depth (fun _ -> random_odd_int ());
    }

let epsilon s = euler /. (float_of_int (1 lsl s.lg_width))
let delta s = 1. /. exp (float_of_int (Array.length s.count))

let update s ~ix ~c =
    Array.iteri (fun i a -> step_matrix s.count i (multiply_shift s.lg_width a ix) c) s.hash_functions

let query s ~ix =
    (* No fusion :-( so this generates an intermediate data structure
     * that is immediately discarded.  We could probably write a
     * minimum_mapi function but meh *)
    minimum (Array.mapi (fun i a -> s.count.(i).(multiply_shift s.lg_width a ix)) s.hash_functions)

(* To implement: *)
(* range_query - needs customized array of sketches *)
(* inner_product_query *)
(* phi_quantiles, heavy_hitters *)

let () =
    let x = make 2.0 0.9 in
    update x 3 4;
    update x 3 2;
    update x 24435 5;
    update x 2323434 1;
    update x 223434 1;
    print_int (query x 234);
    print_string "\n";
    print_float (epsilon x);
    print_string "\n";
    print_float (delta x);
    print_string "\n";
    ()
