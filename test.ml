let rec repeat n thunk =
    if n == 0 then ()
    else (thunk (); repeat (n-1) thunk)

(* RNG initialization and debugging output.  You can modify the
 * generator and seed using GSL_RNG_TYPE and GSL_RNG_SEED.  We
 * reseed OCaml's RNG with this, so our program is entirely
 * deterministic based on these parameters. *)

let () =
    Gsl_error.init ();
    Gsl_rng.env_setup ()

let rng = Gsl_rng.make (Gsl_rng.default ())

let () =
  Printf.printf "\027[34m";
  Printf.printf "gsl rng type=%s seed=%nu\n" (Gsl_rng.name rng) (Gsl_rng.default_seed());
  let seed = Nativeint.to_int (Gsl_rng.get rng) in
  Random.init seed;
  Printf.printf "ocaml rng seed=%d (first gsl rng output)\n" seed;
  Printf.printf "\027[0m";
  Printf.printf "\n"

(** Creates a sketch with parameters [epsilon] and [delta], increments
    a random key (determined by the [kf] thunk) [n] times, and then
    returns a tuple of this sketch and a hash table contanining the
    true counts of the values. *)
let generate ~epsilon ~delta ~n ~kf =
    let s = Cminsketch.make epsilon delta in
    let t = Hashtbl.create (n/8) in (* Worst case needs to grow 3 times *)
    (* It's not necessary to test higher increments, since they
     * are equivalent to calling the function that many times
     * (so we punt to the underling random distribution in kf). *)
    repeat n (fun () -> let k = kf () in
        Cminsketch.update s k 1;
        Hashtbl.replace t k (try Hashtbl.find t k + 1 with Not_found -> 1)
        );
    (s, t)

(** List of distributions we'll test with.  These functions are suitable
    to be passed in as the [kf] parameter in [generate]. *)
let distributions = [
    fun () -> Gsl_randist.poisson rng 10.0;
    ]
let some_distribution = fun () -> Nativeint.to_int (Gsl_rng.get rng)

let () =
    let n = 50000
    and input_epsilon = 0.0001
    and input_delta = 0.8 in
    let (s,t) = generate input_epsilon input_delta n some_distribution in
    let epsilon = Cminsketch.epsilon s
    and delta = Cminsketch.delta s in
    let margin = float_of_int n *. epsilon in
    Printf.printf "epsilon=%.5f delta=%.5f norm=%n, margin=%.0f\n" epsilon delta n margin;
    Printf.printf "width=%d depth=%d size=%d\n" (Cminsketch.width s) (Cminsketch.depth s) (Cminsketch.size s);
    let wrong =
        let wrongref = ref 0 in
        Hashtbl.iter
            (fun k x ->
                let y = Cminsketch.query s k in
                let f = float_of_int (abs (x - y)) in
                if f > margin then wrongref := !wrongref + 1
                    else ()
                )
            t;
        !wrongref in
    let domain_size = Hashtbl.length t in
    Printf.printf "domain=%d\n" domain_size;
    if Cminsketch.size s > domain_size
        then Printf.printf "\027[33mWARNING: Sketch is bigger than hash table\027[0m\n"
        else ();
    let wrong_bound = delta
    and wrong_rate = (float_of_int wrong /. float_of_int domain_size) in
    Printf.printf "wrong: bound=%.5f actual=%.5f\n" wrong_bound wrong_rate;
    Printf.printf "\n";
    (* some statistics about the variance from true value would be cool *)
    if wrong_bound > wrong_rate
        then Printf.printf "\027[1m\027[32mOK\027[0m\n"
        else Printf.printf "\027[1m\027[31mFAIL\027[0m\n"
