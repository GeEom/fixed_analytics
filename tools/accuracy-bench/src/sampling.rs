//! Test point generation.

pub struct SampleStrategy {
    pub grid_points: usize,
    pub random_points: usize,
    pub boundary_points: usize,
    seed: u64,
}

impl SampleStrategy {
    pub fn thorough() -> Self {
        Self {
            grid_points: 5000,
            random_points: 50000,
            boundary_points: 2000,
            seed: 0xDEAD_BEEF_CAFE_BABE,
        }
    }

    pub fn generate(&self, lo: f64, hi: f64) -> Vec<f64> {
        let mut points = Vec::with_capacity(
            self.grid_points + self.random_points + self.boundary_points * 2 + 10,
        );

        for &v in &[0.0, 1.0, -1.0, 0.5, -0.5, 2.0, -2.0] {
            if v >= lo && v <= hi {
                points.push(v);
            }
        }

        for i in 0..self.grid_points {
            let t = i as f64 / (self.grid_points - 1).max(1) as f64;
            points.push(lo + t * (hi - lo));
        }

        let mut rng = self.seed;
        for _ in 0..self.random_points {
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let t = (rng as f64) / (u64::MAX as f64);
            points.push(lo + t * (hi - lo));
        }

        for i in 0..self.boundary_points {
            let t = i as f64 / self.boundary_points.max(1) as f64;
            let exp_t = (1.0 - (-3.0 * t).exp()) / (1.0 - (-3.0_f64).exp());
            let delta = (hi - lo) * 0.1 * (1.0 - exp_t);
            points.push(lo + delta);
            points.push(hi - delta);
        }

        points.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        points.dedup_by(|a, b| (*a - *b).abs() < 1e-15);
        points
    }
}
