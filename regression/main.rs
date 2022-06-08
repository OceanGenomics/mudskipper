use std::process::Command;

fn main() {

        let status = Command::new("bash")
                        .arg("regression/regression_testing.sh")
                        .status()
                        .expect("failed to execute process");
        println!("process finished with: {status}");

        assert!(status.success());
        
}
