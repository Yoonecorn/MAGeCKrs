fn main() {
    env_logger::init();
    if let Err(err) = mageckrs::run() {
        eprintln!("{err}");
        std::process::exit(1);
    }
}
