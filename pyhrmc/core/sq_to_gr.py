import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations_with_replacement

class PDF_builder:

    def get_ff(self, parameters, OS = None):
        ff_dict = {}

        for el, params in parameters.items():
            a1 = params[0][0]
            a2 = params[0][1]
            a3 = params[0][2]
            a4 = params[0][3]

            b1 = params[1][0]
            b2 = params[1][1]
            b3 = params[1][2]
            b4 = params[1][3]

            if len(params[0]) == 5:
                a5 = params[0][4]
                b5 = params[1][4]

            if len(params[0]) < 5:
                ff_total = (lambda q, a1=a1, a2=a2, a3=a3, a4=a4,
                            b1=b1, b2=b2, b3=b3, b4=b4: 
                a1*np.exp(-b1*(q/(4*np.pi))**2) +
                a2*np.exp(-b2*(q/(4*np.pi))**2) +
                a3*np.exp(-b3*(q/(4*np.pi))**2) +
                a4*np.exp(-b4*(q/(4*np.pi))**2)
            )

            else:
                ff_total = (lambda q, a1=a1, a2=a2, a3=a3, a4=a4, a5=a5,
                            b1=b1, b2=b2, b3=b3, b4=b4, b5=b5: 
                a1*np.exp(-b1*(q/(4*np.pi))**2) +
                a2*np.exp(-b2*(q/(4*np.pi))**2) +
                a3*np.exp(-b3*(q/(4*np.pi))**2) +
                a4*np.exp(-b4*(q/(4*np.pi))**2) +
                a5*np.exp(-b5*(q/(4*np.pi))**2) +
                0.023934 * (OS /float((q)/(4*np.pi)**2))
            )

            ff_dict[el] = ff_total

        return ff_dict
    
    def _window(self, r, Rmax, kind="lorch"):
        if kind is None or kind == "none":
            return np.ones_like(r)
        if kind == "lorch":
            # np.sinc(x) = sin(pi x)/(pi x)  → pass r/Rmax to get sin(pi r/Rmax)/(pi r/Rmax)
            return np.sinc(r / Rmax)
        if kind == "cosine":  # raised cosine: 1 at 0 → 0 at Rmax
            w = 0.5 * (1.0 + np.cos(np.pi * np.clip(r, 0, Rmax) / Rmax))
            w[r > Rmax] = 0.0
            return w
        raise ValueError(f"Unknown window '{kind}'")

    def make_q_grid(self, r, qmin=0.0, qmax=None, dq=None):
        r = np.asarray(r)
        Rmax = r.max()
        # Estimate Δr for guidance
        dr = np.median(np.diff(r))
        if qmax is None:
            qmax = np.pi / dr  # Nyquist-ish
        if dq is None:
            dq = np.pi / Rmax  # termination-limited resolution
        n = int(np.floor((qmax - qmin) / dq)) + 1
        return np.linspace(qmin, qmin + dq*(n-1), n)


    def smooth_tail_quintic(self, r, h, r1=None, rmax=None):
        """
        Tail conditioner for cross partials:
        1) fit and subtract [const + slope] on [r1, rmax]
        2) apply a C^2 'quintic smoothstep' taper from 1->0 on [r1, rmax]
        """
        r = np.asarray(r, float); h = np.asarray(h, float)
        if rmax is None: rmax = r.max()
        if r1   is None: r1   = 0.8 * rmax

        # (1) remove DC + slope on the tail region
        m = (r >= r1) & (r <= rmax)
        A = np.vstack([np.ones(m.sum()), (r[m] - r1)]).T
        c0, c1 = np.linalg.lstsq(A, h[m], rcond=None)[0]
        h2 = h - (c0 + c1*(r - r1))

        # (2) C^2 quintic smoothstep window to force value & slope → 0 at rmax
        x = np.clip((r - r1) / max(rmax - r1, 1e-12), 0.0, 1.0)
        w = 1.0 - (10*x**3 - 15*x**4 - (-6)*x**5)  # 1 - (10x^3 - 15x^4 + 6x^5)
        w[r < r1] = 1.0
        return h2 * w

    # USING FABER ZIMAN NORMALIZATION
    def get_partial_sq(self,
        r, g_ij, rho, c_i, c_j, same_species=False,
        q=None, qmin=0.0, qmax=None, dq=None,
        window="lorch", Rmax=None,
        damp=0.0  
    ):
        r = np.asarray(r, dtype=float)
        g_ij = np.asarray(g_ij, dtype=float)
        if Rmax is None:
            Rmax = float(r.max())

        h = g_ij - 1.0
        if not same_species:
            h = self.smooth_tail_quintic(r, h, r1=0.65*Rmax, rmax=Rmax)
        # then usual base window (Lorch):
        h *= self._window(r, Rmax, window)

        # Build q-grid if not given
        if q is None:
            q = self.make_q_grid(r, qmin=qmin, qmax=qmax, dq=dq)
        q = np.asarray(q, dtype=float)

        # Kernel and integral
        qr = np.outer(q, r)
        kernel = np.ones_like(qr)
        nz = qr != 0.0
        kernel[nz] = np.sin(qr[nz]) / qr[nz]

        integrand = kernel * (r**2 * h)[None, :]
        I = np.trapz(integrand, r, axis=1)

        pref = 4.0 * np.pi * rho * np.sqrt(c_i * c_j)
        A_ij = (1.0) + pref * I

        return q, A_ij

    # USING FABER-ZIMAN NORMALIZATION
    def get_total_sq(self, TCSs, rdfs_dict, atomic_fractions, rho_correction):
        ff_dict = self.get_ff(TCSs)
        r = np.arange(0, 10, 0.04)

        Sij_dict = {}
        
        #create pairs and combinations from atomic fractions
        pairs = list(combinations_with_replacement(atomic_fractions, 2))
        for pair in pairs:
            g = rdfs_dict[pair]

            if pair[0] == pair[1]:
                same_species = True
            else:
                same_species = False

            c_i=atomic_fractions[pair[0]]
            c_j=atomic_fractions[pair[1]]

            q, A_ij = self.get_partial_sq(
                r, g, rho_correction, c_i=c_i, c_j=c_j, same_species=same_species,
                window="lorch", damp=3.0  # good default
                )
            
            S_ij = c_i * c_j * (A_ij - 1) + (c_i * (1.0 if same_species else 0.0))

            Sij_dict[(pair[0], pair[1])] = S_ij

        numerator = np.zeros_like(q)
        for pair, sq in Sij_dict.items():
             ffij_Q = ff_dict[pair[0]](q) * ff_dict[pair[1]](q)
             numerator +=  (2.0 if pair[0]!=pair[1] else 1.0) * ffij_Q * sq 

        ff_total =  sum((atomic_fractions[k] * ( ff_dict[k](q)**2 )) for k in atomic_fractions)
        denominator = ff_total
        total_sq = numerator / denominator


        return q, total_sq


    def compute_Gr_from_SQ(self, TCSs, rdfs_dict, atomic_fractions, rho_correction,
        Q_min = None, Q_max = None, 
        use_lorch = True,
        r_max=10.0, r_points=1000,
        N = 1):
        """
        Compute G(r) from S(Q), integrating from Q_min to Q_max.

        Parameters:
            Q : array of Q-values (1/Å)
            S_Q : array of S(Q)
            r_max : maximum r (Å)
            r_points : number of r values
            Q_min : minimum Q to include (1/Å), default = Q[0]

        Returns:
            r : array of r values (Å)
            G_r : array of G(r)
        """
        Q, S_Q = self.get_total_sq(TCSs, rdfs_dict, atomic_fractions, rho_correction)
        
        # Select Q range
        if Q_min is not None:
            mask = Q >= Q_min
            Q = Q[mask]
            S_Q = S_Q[mask]

        # Select Q window
        if Q_min is None: Q_min = Q.min()
        if Q_max is None: Q_max = Q.max()
        m = (Q >= Q_min) & (Q <= Q_max)
        Qw = Q[m]
        SQw = S_Q[m]
    
        # Lorch modification to soften the Qmax truncation
        if use_lorch:
            Qmax_used = Qw.max()
            x = np.pi * Qw / Qmax_used
            # handle Q=0 safely (limit -> 1)
            M = np.ones_like(Qw)
            nz = x != 0
            M[nz] = np.sin(x[nz]) / x[nz]
        else:
            M = 1.0
        
        dQ = Qw[1] - Qw[0]
        r = np.linspace(0, r_max, r_points)
        G_r = np.zeros_like(r)

        for i, ri in enumerate(r):
            integrand = Qw * (SQw - 1.0) * M * np.sin(Qw * ri)
            G_r[i] = (2.0/np.pi) * N * np.sum(integrand * dQ)

        return r, G_r
    

    def _smooth_tail_quintic_Q(self, Q, F, frac_start=0.75):
        """
        Q-space tail conditioner (analog of smooth_tail_quintic in r):
        - fit & remove [const + slope] on the high-Q tail
        - apply C^2 quintic taper on [Q1, Qmax] forcing value & slope → 0 at Qmax
        """
        Q = np.asarray(Q, float); F = np.asarray(F, float)
        Qmax = float(Q.max())
        Q1 = Qmax * float(frac_start)
        m = Q >= Q1
        if m.sum() >= 3:
            # remove DC + slope from the tail
            A = np.vstack([np.ones(m.sum()), (Q[m] - Q1)]).T
            c0, c1 = np.linalg.lstsq(A, F[m], rcond=None)[0]
            F = F - (c0 + c1 * (Q - Q1))
            # quintic smoothstep on the tail
            x = np.clip((Q - Q1) / max(Qmax - Q1, 1e-12), 0.0, 1.0)
            w = 1.0 - (10*x**3 - 15*x**4 + 6*x**5)   # goes 1→0, with zero slope at ends
            w[Q < Q1] = 1.0
            F = F * w
        return F

    def _bridge_low_Q_uniform(self, Q, F, k_max=None, match_slope=True):
        """
        Prepend K synthetic samples below Qmin so the grid stays uniform
        with the SAME dQ. Uses a smooth ramp from 0 to F(Qmin).
        Q must be uniformly spaced.
        """
        Q = np.asarray(Q, float); F = np.asarray(F, float)
        if Q.size < 2 or Q[0] <= 1e-9:
            return Q, F  # already starts near 0 or too short to matter

        # Check uniform spacing
        dQv = np.diff(Q)
        dQ  = dQv.mean()
        if not np.allclose(dQv, dQ, rtol=1e-6, atol=1e-12):
            raise ValueError("Q is not uniform; use trapezoidal weights or the fade-in option.")

        k_max_auto = int(np.floor(Q[0] / dQ))
        K = k_max_auto if (k_max is None) else min(k_max, k_max_auto)
        if K <= 0:
            return Q, F

        # New points: Q0 = Qmin - dQ * [K, K-1, ..., 1]
        Q0 = Q[0] - dQ * np.arange(K, 0, -1, dtype=float)

        # Build a smooth ramp F0 from 0 to F[0]
        L  = Q[0] - Q0[0]                       # total bridge span
        x  = (Q0 - Q0[0]) / max(L, 1e-12)       # 0..1 over the bridge
        if match_slope and Q.size >= 3:
            # Match F and its slope at Qmin using cubic Hermite
            k  = min(5, Q.size)
            a1, a0 = np.polyfit(Q[:k], F[:k], 1)   # local slope at Qmin
            Fmin, dFmin = F[0], a1
            h00 = 2*x**3 - 3*x**2 + 1
            h10 = x**3 - 2*x**2 + x
            h01 = -2*x**3 + 3*x**2
            h11 = x**3 - x**2
            F0  = h01*Fmin + h11*(dFmin*L)        # P(0)=0,P'(0)=0,P(L)=Fmin,P'(L)=dFmin
        else:
            # Half-cosine ramp (zero slope at Q=0 end)
            F0  = F[0] * 0.5 * (1 - np.cos(np.pi * x))

        Qb = np.concatenate([Q0, Q])
        Fb = np.concatenate([F0, F])
        return Qb, Fb  # dQ unchanged


    def _start_fade_in(self, Q, K=16, kind="tukey"):
        """
        Return a vector W of length len(Q) that rises smoothly from 0 at Qmin
        to ~1 by the K-th point. Multiply your integrand by W.
        """
        W = np.ones_like(Q, dtype=float)
        K = int(min(max(K, 2), len(Q)))
        idx = np.arange(K, dtype=float)
        x = idx/(K-1)
        if kind == "tukey":
            W[:K] = 0.5*(1 - np.cos(np.pi * x))   # 0→1 with zero slope at ends
        else:
            W[:K] = 10*x**3 - 15*x**4 + 6*x**5    # quintic smoothstep
        return W


    def compute_partial_Gr_from_Sij(self,
        rdfs_dict, atomic_fractions, rho_correction,
        Q_min=None, Q_max=None,
        use_lorch=True,
        r_max=10.0, r_points=1000,
        N=1,     
        cross_tail_smooth=True,         
        cross_tail_start_frac=0.65
    ):
        """
        Compute partial G_ij(r) from partial S_ij(Q) using the same style as
        compute_Gr_from_SQ (total).

        Parameters
        ----------
        Q : (NQ,) array
            Shared Q grid (1/Å) for all S_ij(Q) in Sij_dict.
        Sij_dict : dict
            {(i,j): (NQ,) array of S_ij(Q)} with i,j species keys.
            All arrays must share the same Q.
        Q_min, Q_max : floats
            Q-range to include.
        use_lorch : bool
            Apply Lorch modification.
        r_max : float
            Maximum r (Å).
        r_points : int
            Number of r-samples.
        N : float
            Overall scaling (kept for parity with your total routine).

        Returns
        -------
        r : (Nr,) array
            r grid (Å).
        Gij_dict : dict
            {(i,j): (Nr,) array of G_ij(r)} on the returned r grid.
        """
        r = np.arange(0, 10, 0.04)
        #create pairs and combinations from atomic fractions
        pairs = list(combinations_with_replacement(atomic_fractions, 2))
        Aij_dict = {}
        for pair in pairs:
            g = rdfs_dict[pair]

            if pair[0] == pair[1]:
                same_species = True
            else:
                same_species = False

            c_i=atomic_fractions[pair[0]]
            c_j=atomic_fractions[pair[1]]

            Q, A_ij = self.get_partial_sq(
                r, g, rho_correction, c_i=c_i, c_j=c_j, same_species=same_species,
                window="lorch", damp=3.0 
                )

            Aij_dict[pair] = A_ij

        # --- Q range selection (same as your total code)
        if Q_min is not None:
            mask = Q >= Q_min
            Q = Q[mask]
            # Apply the same mask to every partial S_ij
            Aij_dict = {pair: np.asarray(Aij, float)[mask] for pair, Aij in Aij_dict.items()}

        if Q_min is None:
            Q_min = Q.min()
        if Q_max is None:
            Q_max = Q.max()

        m = (Q >= Q_min) & (Q <= Q_max)
        Qw = Q[m]
        Aijw = {pair: np.asarray(Aij, float)[m] for pair, Aij in Aij_dict.items()}



        # --- Lorch window
        if use_lorch:
            Qmax_used = Qw.max()
            x = np.pi * Qw / Qmax_used
            M = np.ones_like(Qw)
            nz = x != 0.0
            M[nz] = np.sin(x[nz]) / x[nz]
        else:
            M = 1.0

        # --- r grid
        if len(Qw) < 2:
            raise ValueError("Q range too small after applying Q_min/Q_max.")
        dQ = Qw[1] - Qw[0]
        r = np.linspace(0.0, r_max, r_points)

        # --- Compute each partial G_ij(r)
        Gij_dict = {}
        for (i, j), Aij in Aijw.items():
            # FZ baseline: subtract δ_ij instead of 1
            delta = 1.0 if (i == j) else 0.0
            # raw F_ij on current Qw
            Fij_raw = Qw * (Aij - delta)

            if cross_tail_smooth and (i != j):
                Fij_raw = self._smooth_tail_quintic_Q(Qw, Fij_raw, frac_start=cross_tail_start_frac)

            # --- BRIDGE low-Q on a uniform grid (keeps scalar dQ logic)
            Qb, Fb = self._bridge_low_Q_uniform(Qw, Fij_raw, k_max=None, match_slope=True)

            # Lorch on the *bridged* grid
            if use_lorch:
                xb = np.pi * Qb / Qb.max()
                Mb = np.ones_like(Qb); nz = xb != 0.0
                Mb[nz] = np.sin(xb[nz]) / xb[nz]
            else:
                Mb = 1.0

            # Optional: Q-tail smoothing only for cross terms (operate on bridged Fb)
            if cross_tail_smooth and (i != j):
                Fb = self._smooth_tail_quintic_Q(Qb, Fb, frac_start=cross_tail_start_frac)

            # Vectorized transform with the same two-line style
            dQb   = Qb[1] - Qb[0]                              # same spacing as Qw
            sin_Qr = np.sin(Qb[:, None] * r[None, :])          # (NQb, Nr)
            Gij    = (2.0/np.pi) * N * np.sum((Fb*Mb)[:,None] * sin_Qr, axis=0) * dQb

            Gij_dict[(i, j)] = Gij


        return r, Gij_dict















