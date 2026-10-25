pub type TestResult = Result<(), fixed_analytics::Error>;

pub struct Lcg(pub u64);

impl Lcg {
    pub const fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let x = self.0;
        (x ^ (x >> 29)).wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ (x >> 32)
    }

    pub fn unit(&mut self) -> f64 {
        fixed::types::U0F64::from_bits(self.next_u64() >> 11 << 11).to_num()
    }

    pub fn range(&mut self, lo: f64, hi: f64) -> f64 {
        (hi - lo).mul_add(self.unit(), lo)
    }
}
