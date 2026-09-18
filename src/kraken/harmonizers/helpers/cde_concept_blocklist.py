"""
Curated blocklist of NCIt/UMLS concepts too generic/grammatical to be meaningful edge targets
(numbers, prepositions, comparators, degree & temporal glue, response-glue like Unknown/Other,
and phrase fragments). Covers both the concept fields (objectClass/property/dataElementConcept ->
assesses edges) and the constituents of compound permissible-value codes (-> related_to edges).

Derivation: each candidate concept was looked up in the integrated kraken node set; those typed
ONLY as biolink:NamedThing were surfaced as junk candidates, and a human kept the substantive
ones (see RESCUE_IDS / RESCUE_ANSWER_NAMES in scripts/gen_cde_concept_blocklist.py). All names below
are pulled verbatim from the kraken node file. Regenerate/review when the CDE export changes.

The harmonizer skips both the edge AND the stub node for any concept in this set.
"""

CONCEPT_BLOCKLIST: frozenset[str] = frozenset(
    {
        "UMLS:C3166761",  #
        "UMLS:C3177045",  #
        "UMLS:C3176995",  #
        "UMLS:C1317484",  #
        "UMLS:C3166424",  #
        "UMLS:C3177086",  #
        "UMLS:C3483359",  #
        "UMLS:C3177042",  #
        "UMLS:C3177061",  #
        "NCIT:C78209",  # ability
        "NCIT:C14140",  # Acute
        "NCIT:C25382",  # Administered
        "NCIT:C50400",  # Age Unit
        "NCIT:C25415",  # Alternative
        "NCIT:C48928",  # And/Or
        "NCIT:C25281",  # associated
        "NCIT:C25427",  # At
        "NCIT:C25213",  # Baseline
        "NCIT:C68615",  # Birth Date
        "NCIT:C85883",  # Burden
        "NCIT:C49152",  # Case
        "NCIT:C25446",  # change
        "NCIT:C14141",  # chronic
        "NCIT:C25161",  # classification
        "NCIT:C37926",  # Cognitive
        "NCIT:C79892",  # Concern
        "NCIT:C25730",  # Concomitant
        "NCIT:C25456",  # concurrent
        "NCIT:C25457",  # Condition
        "NCIT:C53279",  # Continue
        "NCIT:C61299",  # Control
        "NCIT:C69088",  # cost
        "NCIT:C25467",  # Critical
        "NCIT:C25471",  # Current
        "NCIT:C142500",  # Data Validation
        "NCIT:C25164",  # date
        "NCIT:C37939",  # Date and Time
        "NCIT:C25301",  # Day
        "NCIT:C25476",  # Delay
        "NCIT:C41205",  # disposition
        "NCIT:C62289",  # Domain
        "NCIT:C61589",  # Early
        "NCIT:C75539",  # Episode
        "NCIT:C16390",  # Etiology
        "NCIT:C64651",  # Every
        "NCIT:C73992",  # excessive
        "NCIT:C102844",  # Extreme
        "NCIT:C25745",  # Failure
        "NCIT:C25509",  # first
        "NCIT:C25511",  # Focus
        "NCIT:C66835",  # Four
        "NCIT:C64649",  # Frequently
        "NCIT:C25516",  # From
        "NCIT:C25517",  # full
        "NCIT:C68846",  # Global
        "NCIT:C61421",  # Greater
        "NCIT:C61584",  # greater than
        "NCIT:C25178",  # Health
        "NCIT:C25227",  # High
        "NCIT:C25529",  # hour
        "NCIT:C25737",  # Identification
        "NCIT:C16726",  # Incidence
        "NCIT:C166400",  # Inclusive
        "NCIT:C25180",  # Indicator
        "NCIT:C121660",  # Influence
        "NCIT:C48364",  # Instance
        "NCIT:C45255",  # integer
        "NCIT:C72886",  # interference
        "NCIT:C14159",  # Invasive
        "NCIT:C25549",  # Isolation
        "NCIT:C25551",  # Last
        "NCIT:C61585",  # Less Than
        "NCIT:C25554",  # Level
        "NCIT:C54722",  # low
        "NCIT:C16833",  # Medicine
        "NCIT:C25569",  # middle
        "NCIT:C29846",  # month
        "NCIT:C128667",  # More
        "NCIT:C64934",  # Morning
        "NCIT:C107222",  # Nearly
        "NCIT:C41204",  # need
        "NCIT:C25594",  # Negation
        "NCIT:C25586",  # New
        "NCIT:C48660",  # Not Applicable
        "NCIT:C126101",  # Not Available
        "NCIT:C49484",  # Not Done
        "NCIT:C43234",  # Not Reported
        "NCIT:C25337",  # number
        "NCIT:C52349",  # Obtain
        "NCIT:C83477",  # Occupancy
        "NCIT:C127786",  # Occurrence Indicator
        "NCIT:C66832",  # One
        "NCIT:C25279",  # Onset
        "NCIT:C37998",  # Or
        "NCIT:C25311",  # oral
        "NCIT:C17649",  # other
        "NCIT:C151896",  # Other Than
        "NCIT:C20200",  # Outcome
        "NCIT:C25605",  # Overall
        "NCIT:C25609",  # Past
        "NCIT:C95402",  # Past Week
        "NCIT:C180612",  # pay
        "NCIT:C65039",  # Per
        "NCIT:C38000",  # Performed
        "NCIT:C25233",  # peripheral
        "NCIT:C43623",  # Persistent
        "NCIT:C25618",  # physical
        "NCIT:C25319",  # Place
        "NCIT:C53344",  # Population-Based Statistic
        "NCIT:C38008",  # Post
        "NCIT:C25629",  # Prior
        "NCIT:C54104",  # Private
        "NCIT:C94316",  # Psychological
        "NCIT:C13304",  # Pulmonary
        "NCIT:C49143",  # Qualitative Evaluation
        "NCIT:C25639",  # Receive
        "NCIT:C25197",  # recommendation
        "NCIT:C25648",  # Relationship
        "NCIT:C25652",  # Requirement
        "NCIT:C17102",  # Risk
        "NCIT:C25666",  # second
        "NCIT:C25200",  # self
        "NCIT:C188309",  # Self-Reported Information
        "NCIT:C107221",  # several
        "NCIT:C25676",  # Severity
        "NCIT:C66837",  # Six
        "NCIT:C62650",  # social
        "NCIT:C65099",  # Some
        "NCIT:C38024",  # Specified
        "NCIT:C25685",  # Specify
        "NCIT:C25688",  # Status
        "NCIT:C25515",  # Temporal Frequency
        "NCIT:C47891",  # Test
        "NCIT:C25277",  # Therapeutic
        "NCIT:C66834",  # Three
        "NCIT:C65107",  # To
        "NCIT:C25304",  # Total
        "NCIT:C66833",  # Two
        "NCIT:C25284",  # Type
        "NCIT:C28012",  # Unilateral
        "NCIT:C25709",  # Unit of Measure
        "NCIT:C17998",  # Unknown
        "NCIT:C64921",  # Upon Awakening
        "NCIT:C95018",  # Use
        "NCIT:C102843",  # Usual
        "NCIT:C94301",  # Vasoactive
        "NCIT:C45513",  # Verification
        "NCIT:C27985",  # Viral
        "NCIT:C29844",  # week
        "NCIT:C25718",  # Without
        "NCIT:C29848",  # year
    }
)
