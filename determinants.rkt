#lang racket

;; A Row is a (List-of Number)
;; Examples:
(module+ test
  (require rackunit)
  (check-equal? (make-row '(1 1 1)) '(1 1 1)))

(define (make-row ls)
  (cond [(empty? ls) (error 'make-row "a row cannot be empty")]
        [else ls]))

;; A Matrix is either
;;  - Empty OR
;;  - (make-matrix (List-of Row))

;; Examples:
;; empty
(module+ test
  (check-equal? (make-matrix) empty)
  (check-equal? (make-matrix '((1 0) (0 1))) '((1 0) (0 1)))
  (check-equal? (make-matrix '((1 0 0 -2) (0 1 0 -3) (0 0 1 -4)))
                '((1 0 0 -2) (0 1 0 -3) (0 0 1 -4)))
  ;; (check-exn (make-matrix '((1 1 1) (1 1)))
  ;;            "error: rows of different lengths")
  )

(define (make-matrix [ls-of-ls empty])
  (local [(define (same-numbers? ls-of-ls)
            (cond [(or (= 0 (length ls-of-ls)) (= 1 (length ls-of-ls))) #t]
                  [else (if (not (= (first ls-of-ls) (second ls-of-ls)))
                            #f
                            (same-numbers? (rest ls-of-ls)))]))]
         (cond [(empty? ls-of-ls) empty]
               [(not (same-numbers? (map length ls-of-ls)))
                (error 'make-matrix "rows of matrix don't have the same lengths")]
               [else ls-of-ls])))

;; Matrix -> (List-of Number) OR Symbol

;; Takes a matrix representing a system of equations Ax=b and produces
;; the solution set for the system.  The vector b must be provided in
;; the matrix (as the last column.)  If the system has a non-unique
;; solution or no solution, I merely produce

;;                'no-solution-or-non-unique-solution

;; because I don't know how to do better.

;; Examples:
(module+ test
  (check-equal? (solve '((1 0 -2) (0 1 -3))) '(2 3))
  (check-equal? (solve '((1 0 0 2) (0 1 0 3) (0 0 1 4))) '(-2 -3 -4))
  (check-equal? (solve '((0 0 0 0) (0 1 0 3) (0 0 1 4))) 'no-solution-or-non-unique-solution)
  (check-equal? (solve '((1 0 -2) (0 x -3))) 'non-numeric-matrix))

(define (solve m)
  (if (not (matrix-numbers-only? m))
      'non-numeric-matrix
      (let* ([rows (matrix-rows m)]
             [cols (matrix-cols m)]
             [Dc (determinant (matrix-erase-column m cols))])
        (if (= Dc 0)
            'no-solution-or-non-unique-solution
            (map (lambda (i)
                   (solve-for m i)) (build-list rows (lambda (t) (add1 t))))))))

;; A VariableIndex is an Index, but with a restriction.  For a matrix
;; (make-matrix '((1 0 -2) (0 1 -3))), while Index can be any of 1, 2
;; and 3, VariableIndex can only be 1 or 2.  In other words,
;; VariableIndex ranges on the unknowns of the system, not on the
;; columns of the matrix.

;; Matrix VariableIndex -> Number

;; Takes a matrix m, an index x specifying a variable to be solved for
;; and produces the value of that variable, if the system has a
;; solution.  We assume the system always has a unique solution.

;; Examples:
(module+ test
  (check-equal? (solve-for '((1 0 -2) (0 1 -3)) 2) 3)
  (check-equal? (solve-for '((1 0 -2) (0 1 -3)) 1) 2)
  (check-equal? (solve-for '((1 0 0 -2) (0 1 0 -3) (0 0 1 -4)) 1) 2)
  (check-equal? (solve-for '((1 0 0 -2) (0 1 0 -3) (0 0 1 -4)) 2) 3)
  (check-equal? (solve-for '((1 0 0 -2) (0 1 0 -3) (0 0 1 -4)) 3) 4))

;; The names Dx and Dc follow up on the notation by Lancelot Thomas
;; Hogben, exposed in ``Mathematics for the Million'', W. W. Norton,
;; 1968.  See the chapter ``The Algebra of the Chessboard.''

;; Essentially, Dx is the determinant relative to the x-unknown.  It
;; is the determinant obtained by ignoring the x-column and
;; considering only the other columns.  The Dc determinant is the
;; determinants of the [c]onstants; that is, the determinant obtained
;; by ignoring the constants.  (Ignoring the b-vector of Ax=b.)

;; Key equation here is:
;; x/Dx = -y/Dy = z/Dz = -1/Dc

;; It means -x = Dx/Dc
;;           y = Dy/Dc
;;          -z = Dz/Dc

;; That's for a 3x4 matrix. If it were a 4x5 matrix, then the key
;; equation would be -x/Dx = y/Dy = -z/Dz = w/Dw = -1/Dc

;; It means x = Dx/Cc
;;         -y = Dy/Dc
;;          z = Dz/Dc
;;         -w = Dw/Dc

;; In other words, an odd number of equations means we begin negative
;; and alternate signs, while even number means we begin positive and
;; alternative signs.

(define (solve-for m i)
  (let* ([rows (matrix-rows m)]
         [cols (matrix-cols m)]
         [Dx (determinant (matrix-erase-column m i))]
         [Dc (determinant (matrix-erase-column m cols))]
         [sign (lambda (i) (expt -1 (if (odd? rows) i (add1 i))))])
    (* (/ Dx Dc) (sign i))))

;; Matrix -> Number

;; Takes a matrix m and produces the determinant of m.

;; Examples:
(module+ test
  (check-equal? (determinant empty) 0)
  (check-equal? (determinant '((9))) 9)
  (check-equal? (determinant '((1 1))) 'non-square)
  (check-equal? (determinant '((1 0) (0 1))) 1)
  (check-equal? (determinant '((1 1) (1 1))) 0)
  (check-equal? (determinant '((1 2 3) (3 5 7) (7 9 9))) 2)
  (check-equal? (determinant '((1 2 3) (4 5 6) (7 8 9))) 0))

(define (determinant m)
  (local
   [(define (determinant-2-by-2 m)
      (let ([a (matrix-ref m 1 1)]
            [b (matrix-ref m 1 2)]
            [c (matrix-ref m 2 1)]
            [d (matrix-ref m 2 2)])
        (- (* a d) (* b c))))
    (define (determinant-any-dimension m)
      (cond [(is-matrix-1-by-1? m) (matrix-ref m 1 1)]
            [(is-matrix-2-by-2? m) (determinant-2-by-2 m)]
            [else
             (for/fold
                 ([acc-sum 0])
                 ([j (in-naturals 1)]
                  [e (in-list (matrix-row-ref m 1))])
               (let ([s (if (even? j) -1 1)]) ; sign of co-factor
                 (+ acc-sum
                    (* s e (determinant-any-dimension (minor m 1 j))))))]))]
   (cond [(empty? m) 0]
         [(not (is-matrix-square? m)) 'non-square]
         [else (determinant-any-dimension m)])))

;; Matrix Index Index -> Matrix

;; Produces the corresponding minor for the co-factor expansion of a
;; determinant.  More precisely, it takes a matrix, a RowIndex, a
;; ColumnIndex, producing a matrix without the row and column
;; specified.

;; Examples:
(module+ test
  (check-equal? (minor '((1 2 3) (3 5 7) (7 9 9)) 1 1) '((5 7) (9 9)))
  (check-equal? (minor '((1 2 3) (3 5 7) (7 9 9)) 1 2) '((3 7) (7 9))))

(define (minor m r c)
  (matrix-erase-row (matrix-erase-column m c) r))

;; Matrix Index -> Matrix

;; It takes a matrix and a column-index, producing a matrix from the
;; former without the specified column.

;; Examples:
(module+ test
  (check-equal?
   (matrix-erase-column '((1 2 3) (3 5 7) (7 9 9)) 1)
   '((2 3) (5 7) (9 9))))

(define (matrix-erase-column m c)
  (cols->rows (untag (filter (lambda (pair) (not (= (car pair) c)))
                             (tag (rows->cols m))))))

;; Matrix Index -> Matrix

;; It takes a matrix and a row-index, producing a matrix from the
;; former without the specified row.

;; Examples:
(module+ test
  (check-equal?
   (matrix-erase-row '((1 2 3) (3 5 7) (7 9 9)) 2)
   '((1 2 3) (7 9 9))))

(define (matrix-erase-row m r)
  (untag (filter (lambda (pair) (not (= (car pair) r))) (tag m))))

;; (List-of (Pair X Y)) -> (List-of Y)

;; Takes a list of pairs, which are composed of types X and Y,
;; producing a list of Y only.  That it, produce a list of the Ys in
;; the pairs, getting rid of the Xs.

;; Examples:
(module+ tets
  (untag '((1 (x x x)) (2 (y y y)))) '((x x x) (y y y)))

(define (untag ls-of-pairs)
  (for/list ([pair ls-of-pairs])
    (cadr pair)))

;; (List-of Y) -> (List-of (Pair Natural Y))

;; Takes a list of Y and produces a list of pairs.  Each pair is
;; composed of a natural number and Y.  That is, index the list out
;; with natural numbers.

;; Examples:
(module+ test
  (check-equal? (tag '((1 2 3) (3 5 7))) '((1 (1 2 3)) (2 (3 5 7))))
  (check-equal? (untag (tag '((1 2 3) (3 5 7)))) '((1 2 3) (3 5 7))))

(define (tag ls-of-y)
  (for/list ([e ls-of-y]
             [i (in-naturals 1)])
    (list i e)))

;; Matrix Natural -> Boolean

;; Takes a matrix m, a natural dim, the dimension desired, producing
;; TRUE whether the matrix is square and is of the desired dimension.

;; Examples:
(module+ test
 (check-equal? (is-matrix-dim-by-dim? '((1 1) (1 1)) 2) #t)
 (check-equal? (is-matrix-dim-by-dim? '((1 1 1)) 2) #f))

(define (is-matrix-dim-by-dim? m dim)
  (= dim (matrix-rows m) (matrix-cols m)))

(define (is-matrix-2-by-2? m)
  (is-matrix-dim-by-dim? m 2))

(define (is-matrix-1-by-1? m)
  (is-matrix-dim-by-dim? m 1))

;; Matrix -> Boolean

;; Takes a matrix, producing TRUE is the matrix is square.

;; Examples:
(module+ test
  (check-equal? (is-matrix-square? '((1 0) (0 1))) #t)
  (check-equal? (is-matrix-square? '()) #t)
  (check-equal? (is-matrix-square? '((9))) #t)
  (check-equal? (is-matrix-square? '((1 1 1) (1 1 1) (1 1 1))) #t)
  (check-equal? (is-matrix-square? '((1 1))) #f))
  
(define (is-matrix-square? m)
  (= (matrix-rows m) (matrix-cols m)))

;; Matrix -> Matrix

;; Takes a row-oriented matrix m and produces a column-oriented
;; matrix.

;; Examples:
(module+ test
  (check-equal? (rows->cols '((1 2 3) (4 5 6))) '((1 4) (2 5) (3 6)))
  (check-equal? (rows->cols '((1 2 3) (4 5 6) (3 3 3))) '((1 4 3) (2 5 3) (3 6 3)))
  (check-equal? (rows->cols empty) empty))

(define (rows->cols m)
  (for/list ([j (in-range 1 (add1 (matrix-cols m)))])
    (for/list ([i (in-range 1 (add1 (matrix-rows m)))])
      (matrix-ref m i j))))

;; Matrix -> Matrix

;; Takes a column-oriented matrix m and produces a row-oriented
;; matrix.  As it turns out, this function is identical to the above.

;; Examples:
(module+ test
  (check-equal? (cols->rows (rows->cols '((1 2 3) (4 5 6)))) '((1 2 3) (4 5 6))))

(define cols->rows rows->cols)

;; Matrix -> Boolean

;; Takes a matrix m, producing TRUE if the matrix is a numeric one.

;; Examples:
(module+ test
  (check-equal? (matrix-numbers-only? (make-matrix empty)) #t)
  (check-equal? (matrix-numbers-only? (make-matrix '((1 1) (1 3)))) #t)
  (check-equal? (matrix-numbers-only? (make-matrix '((1 x)))) #f))

(define (matrix-numbers-only? m)
  (local [(define (true? x) (not (false? x)))]
         (cond [(empty? m) #t]
               [else (andmap true? (for/list ([row m])
                                     (andmap number? row)))])))

;; Matrix Index -> Row

;; Takes a matrix m, a row-index i, procuding the i-th row of m

;; Examples:
(module+ test
  (check-equal? (matrix-row-ref '((1 2 3) (3 5 7) (7 9 9)) 2) '(3 5 7)))

(define (matrix-row-ref m r)
  (list-ref m (sub1 r)))

;; Matrix Index Index -> Number

;; Takes a matrix and two indices, producing the Matrix(i,j) element

;; Examples:
(module+ test
  (check-equal? (matrix-ref '((1 2 3) (5 7 8) (4 6 3)) 1 1) 1)
  (check-equal? (matrix-ref '((1 2 3) (5 7 8) (4 6 3)) 2 2) 7))

(define (matrix-ref matrix i j)
  (list-ref (list-ref matrix (sub1 i)) (sub1 j)))

;; Matrix -> Natural

;; Takes a matrix producing its number of rows

;; Examples:
(module+ test
  (check-equal? (matrix-rows '((1 0 2) (0 1 3))) 2)
  (check-equal? (matrix-rows empty) 0))

(define (matrix-rows m)
  (length m))

;; Matrix -> Natural

;; Takes a matrix producing its number of columns

;; Examples:
(module+ test
  (check-equal? (matrix-cols '((1 0 2) (0 1 3))) 3))

(define (matrix-cols m)
  (if (empty? m) 0 (length (first m))))

(provide (all-defined-out))

(module+ main
  (printf "~a.\n" "Hello world"))
