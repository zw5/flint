"""Numerical regression checks. Run: python -m unittest discover -s tests."""

import csv
import io
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import cross_validate_readouts as cv
from verify import benjamini_hochberg, partial_spearman


def csv_text(rows: list[dict]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


class CalibrationTests(unittest.TestCase):
    def test_held_out_values_do_not_enter_the_fit(self):
        predictors = np.array([[0.0], [1.0], [2.0], [3.0], [1e9]])
        outcome = np.array([2.0, 5.0, 8.0, 11.0, -1e12])
        train, test = np.arange(4), np.array([4])
        prediction = cv.fit_predict(predictors, outcome, train, test)
        np.testing.assert_allclose(prediction, [3e9 + 2], rtol=1e-14)
        outcome[4] = 1e12
        np.testing.assert_array_equal(
            prediction, cv.fit_predict(predictors, outcome, train, test)
        )

    def test_no_predictors_uses_training_mean(self):
        outcome = np.array([2.0, 4.0, 6.0, 1000.0])
        prediction = cv.fit_predict(
            np.empty((4, 0)), outcome, np.arange(3), np.array([3])
        )
        np.testing.assert_array_equal(prediction, [4.0])

    def test_each_plan_holds_every_participant_once(self):
        plans = cv.make_fold_plans(111)
        self.assertEqual(len(plans), 21)
        for plan in plans:
            held = [index for fold in plan.test_folds for index in fold]
            self.assertEqual(sorted(held), list(range(111)))
        repeated = cv.make_fold_plans(111)
        self.assertEqual(plans, repeated)


class StatisticalTests(unittest.TestCase):
    def test_fdr_restores_original_order(self):
        np.testing.assert_allclose(
            benjamini_hochberg([0.01, 0.04, 0.03, 0.002]), [0.02, 0.04, 0.04, 0.008]
        )

    def test_partial_ranks_drop_incomplete_participants(self):
        score = np.array([1.0, 2.0, 3.0, 4.0, np.nan, 8.0])
        outcome = np.array([2.0, 1.0, 4.0, 3.0, 5.0, 6.0])
        control = np.array([0.0, 1.0, 0.0, 1.0, 0.0, np.nan])
        actual = partial_spearman(score, outcome, [control])
        expected = partial_spearman(score[:4], outcome[:4], [control[:4]])
        self.assertEqual(actual, expected)
        self.assertEqual(actual[2], 4)


class PublishedResultTests(unittest.TestCase):
    def test_all_folds_predictions_and_metrics_match_published_tables(self):
        participant_ids, readouts = cv.load_readouts()
        plans = cv.make_fold_plans(len(participant_ids))
        predictions, metrics, checks = cv.evaluate_readouts(
            participant_ids, readouts, plans
        )
        summaries = cv.summarize_runs(metrics, readouts, len(participant_ids))
        self.assertEqual(checks, 12)
        tables = {
            "fold_assignments.csv": cv.fold_records(plans, participant_ids),
            "predictions.csv": predictions,
            "run_metrics.csv": metrics,
            "summary.csv": summaries,
        }
        for filename, rows in tables.items():
            with self.subTest(filename=filename):
                expected = (ROOT / "results/cross_validation" / filename).read_bytes()
                self.assertEqual(csv_text(rows).encode(), expected)


if __name__ == "__main__":
    unittest.main()
