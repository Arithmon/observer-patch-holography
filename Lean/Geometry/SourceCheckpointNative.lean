import Geometry.SourceCheckpointMeaning
import Geometry.SourceCheckpointPolicy
import Geometry.SourceCheckpointEntropy

/-!
# Composition of native observability, finite feasibility and checkpoint selection

The transition state is a response checkpoint over all admissible
preparations, not a realized payload. Consequently the same transition law
applies to every source intervention. Topology supplies a feasible finite
deadline and the backward recursion supplies a terminating observation law.
Selecting the actual public record interface and physical deadline remains
a source obligation.
-/

set_option autoImplicit false

namespace OPH.SourceCheckpointNative
noncomputable section
open OPH.SourceTemporalGuard OPH.SourceCheckpointMeaning OPH.SourceCheckpointPolicy

variable {X R P A : Type*} [DecidableEq P] [Fintype A]

def step (root : P) (seam : A → P × P) (c : Checkpoint X P) (a : A) : Checkpoint X P :=
  advance root c (seam a)

def terminal (meaning : X → R) (c : Checkpoint X P) : ℝ := by
  classical
  exact if Published meaning c then 1 else 0

theorem terminal_nonneg (meaning : X → R) (c : Checkpoint X P) : 0 ≤ terminal meaning c := by
  classical
  unfold terminal
  split <;> norm_num

theorem terminal_positive_iff (meaning : X → R) (c : Checkpoint X P) :
    0 < terminal meaning c ↔ Published meaning c := by
  classical
  by_cases h : Published meaning c <;> simp [terminal,h]

theorem finish_run (root : P) (seam : A → P × P) (c : Checkpoint X P) (as : List A) :
    finish (step root seam) c as = run root (as.map seam) c := by
  induction as generalizing c with
  | nil => rfl
  | cons a as ih => exact ih _

theorem native_positive_iff (meaning : X → R) (root : P) (seam : A → P × P)
    (weight : Checkpoint X P → A → ℝ) (hw : ∀ c a, 0 < weight c a)
    (c : Checkpoint X P) (as : List A)
    (hm : 0 < mass (step root seam) weight (terminal meaning) as.length c) :
    0 < pathLaw (step root seam) weight (terminal meaning) c as ↔
      Published meaning (run root (as.map seam) c) := by
  rw [pathLaw_positive_iff _ _ _ hw (terminal_nonneg meaning) c as hm,
    terminal_positive_iff,finish_run]

/-- Each prefix with positive continuation mass retains exactly the public
meaning needed for continuation. Irrelevant source distinctions may vanish. -/
theorem positive_mass_viable (meaning : X → R) (root : P) (seam : A → P × P)
    (weight : Checkpoint X P → A → ℝ) (hw : ∀ c a, 0 < weight c a)
    (c : Checkpoint X P) (n : ℕ)
    (hm : 0 < mass (step root seam) weight (terminal meaning) n c) : Viable meaning c := by
  obtain ⟨as, _, ha⟩ := (mass_positive_iff _ _ _ hw (terminal_nonneg meaning) n c).mp hm
  rw [finish_run,terminal_positive_iff] at ha
  exact publication_requires_viability meaning root c (as.map seam) ha

theorem encode_supported (G : SimpleGraph P) (root : P) (seam : A → P × P)
    (cover : ∀ e : P × P, G.Adj e.1 e.2 →
      ∃ a, ∀ c : Checkpoint X P, advance root c (seam a) = advance root c e)
    (es : List (P × P)) (hs : SourceTemporalConnected.Supported G es) :
    ∃ as : List A, as.length = es.length ∧
      ∀ c : Checkpoint X P, run root (as.map seam) c = run root es c := by
  induction es with
  | nil => exact ⟨[],rfl,fun _ => rfl⟩
  | cons e es ih =>
    obtain ⟨a,ha⟩ := cover e (hs e (by simp))
    obtain ⟨as,hl,he⟩ := ih (fun e he => hs e (by simp [he]))
    refine ⟨a::as,by simp [hl],?_⟩
    intro c
    simp only [List.map_cons,run,ha,he]

/-- Connected topology removes a feasibility assumption: a viable sampled
checkpoint has a positive finite partition function at a bounded horizon. -/
theorem connected_positive_deadline [Fintype P] (G : SimpleGraph P) (hc : G.Connected)
    (meaning : X → R) (root : P) (seam : A → P × P)
    (cover : ∀ e : P × P, G.Adj e.1 e.2 →
      ∃ a, ∀ c : Checkpoint X P, advance root c (seam a) = advance root c e)
    (weight : Checkpoint X P → A → ℝ) (hw : ∀ c a, 0 < weight c a)
    (c : Checkpoint X P) (hv : Viable meaning c) (hr : HasRoot root c) :
    ∃ n, 2*n + Fintype.card P ≤ Fintype.card P ^ 2 ∧
      0 < mass (step root seam) weight (terminal meaning) n c := by
  obtain ⟨es,hs,hb,hp⟩ := (connected_viability_iff G hc meaning root c hr).mp hv
  obtain ⟨as,hl,he⟩ := encode_supported G root seam cover es hs
  refine ⟨as.length,by simpa [hl] using hb,?_⟩
  apply positive_path_mass _ _ _ hw (terminal_nonneg meaning) c as
  rw [finish_run,he,terminal_positive_iff]
  exact hp

end
end OPH.SourceCheckpointNative
