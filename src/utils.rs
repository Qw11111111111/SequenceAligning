use std::char;

pub fn vec_u8_to_str(data: &[u8]) -> String {
    data.iter().map(|&i| i as char).collect::<String>()
}
