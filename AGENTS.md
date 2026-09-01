# Instructions for this R package

Whenever you create or modify a function in `R/`:

1. Create or update a corresponding test in `tests/testthat/`.
2. Use the `testthat` framework.
3. Test:
   - normal inputs,
   - edge cases,
   - invalid inputs and expected errors.
4. Run the relevant tests.
5. Run the complete test suite with:

   ```sh
   Rscript -e "devtools::test()"
   ```

6. Before finishing, run:

   ```sh
   R CMD check .
   ```

7. Do not say the task is complete if tests fail.
8. Report which tests were added and whether they passed.
