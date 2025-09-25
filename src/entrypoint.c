// We need to forward routine registration from C to Rust
// to avoid the linker removing the static library.

void R_init_sovlephrust_extendr(void *dll);

void R_init_sovlephrust(void *dll) {
    R_init_sovlephrust_extendr(dll);
}
