use std::env;

mod annotations;

fn main() {
    let src = env::args().nth(1).expect("missing src");
    annotations::read(src);
}
