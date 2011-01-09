(** This module implements the count-min sketch, a sublinear
    space, probabilistic data structure
    invented by Graham Cormode and S Muthukrishnan, described in
    "An Improved Data Stream Summary: The Count-Min Sketch and its
    Applications."  It is well suited for summarizing data streams and
    finding quantiles/frequent items.  It has also found novel
    uses: see "Popularity is Everything" by Schechter, Herley
    and Mitzenmacher for an approach that uses the count min sketch
    to protect passwords from statistical guessing attacks.

    This implementation is presently incomplete: it still needs
    more of the original functions described in the original paper,
    functionality and statistical tests and a serialization
    format.  Future directions include incorporating the saturation
    protection mechanisms described in "The Eternal Sunshine of the
    Sketch Data Structure." *)

type sketch

(** Multiply shift, a weak universal hash family that this implemenation
    uses to back its sketches (a more conventional choice that was not
    used in this implementation is Carter
    and Wegman's universal hash family).  It was presented in
    "A Reliable Randomized Algorithm for the Closest-Pair Problem"
    and has the property that for all ints [x] and [y] such that
    [x != y], and for [h] chosen uniformly at random from [H =
    {multiply_shift m a | a is odd}], which maps ints to ints of range
    [0] to [2^m - 1],

    {[ Pr[h(x) = h(y)] <= 1/m ]}

    A notable quirk about this hash function is that the size of its
    output range must be a power of two.
    *)
val multiply_shift : m:int -> a:int -> x:int -> int

(** Create a count-min sketch for which the error in answering
    a query is within a factor of [epsilon], and gets it totally wrong
    only with probability [delta].  You get more accurate results for
    small epsilon and small delta, but use less memory for larger
    epsilon and larger delta.

    More detailed bounds regarding [epsilon] and [delta] can be found
    in the relevant estimation functions.

    Side effects: Uses the random number generator to generate
    the hash functions.

    @raise Invalid_argument if [epsilon <= 0.0], [delta >= 1.0] or
           [delta <= 0.0]. *)
val make : epsilon:float -> delta:float -> sketch

(** Creates a count-min sketch, manually specifying the dimensions of
    the backing matrix.  The relationship between [depth]/[width] and
    [epsilon]/[delta] is as follows: [width = ceil(e / epsilon)] and
    [depth = ceil(ln(1/delta))] (where [e] denotes Euler's constant).
    You can determine the effective [epsilon] and [delta] using the
    respective functions.

    @raise Invalid_argument if [depth <= 0] or [width <= 0]. *)
val make_int : depth:int -> width:int -> sketch

(** Returns the true error factor for a sketch. *)
val epsilon : sketch -> float

(** Returns the true error probability for a sketch. *)
val delta : sketch -> float

(** Returns the width of the sketch. *)
val width : sketch -> int

(** Returns the depth of the sketch. *)
val depth : sketch -> int

(** Returns the size of the sketch.  This should be smaller than your
    input space, otherwise you're not getting any benefit! *)
val size : sketch -> int

(** [O(log(1/delta))]. Updates a sketch adding [c] to the field [ix]. *)
val update : sketch -> ix:int -> c:int -> unit

(** [O(log(1/delta))]. Estimates the count of the field [ix] if all the
    actual counts are non-negative.  Use [nquery] if actual counts may b
    negative.

    This estimate is never less than the true value and, with
    probability of at least [1 - delta], the overestimation is no
    greater than [epsilon * |a|1], where [|a|1] denotes the L1 (taxicab)
    norm on the actual vector [a] (i.e. the sum of all updates done to
    all keys in the sketch.) *)
val query : sketch -> ix:int -> int

(** [~O(log(1/delta))]. Estimates the count of the
    field [ix].  This works in both negative and non-negative actual
    counts, but is not as accurate as [query] in the non-negative case.

    With probability of at least [1 - delta^(1/4)], the estimate falls
    within [3  epsilon * |a|1] of the true value.

    The true asymptotic time complexity of this particular
    implementation is actually [O(log(1/delta) * log(log(1/delta))]
    due to the use of an [O(n * log(n)] median algorithm. *)
val nquery : sketch -> ix:int -> int

(** [O(1/epsilon) * log(1/delta)]. Estimates the dot product of the two
    sketches, that is, [sum_i(a.(i) *  b.(i))], when all actual counts
    in both sketches are non-negative.

    This estimate is never less than the true value and, with
    probability of at least [1 - delta], the overestimation is no
    greater than [epsilon * |a|1 * |b|1].  [query] is in fact
    a specialiation of this operation in the case that [b] is
    a unit vector. *)
val dot_product_query : sketch -> sketch -> int

(** Estimates the dot product the two sketches for negative and
    non-negative cases.  As with [nquery], this suffers an
    accuracy hit, and furthermore, there are no known practical
    uses for the dot product in the case of negative vectors.
    Included for completeness. *)
val dot_product_nquery : sketch -> sketch -> int
