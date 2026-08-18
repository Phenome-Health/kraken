# Changelog

## [2.1.2](https://github.com/Phenome-Health/kraken/compare/v2.1.1...v2.1.2) (2026-08-18)


### Documentation

* clarify python version compatibility ([b5d0c83](https://github.com/Phenome-Health/kraken/commit/b5d0c83a326a58c65f53b5049e54378a65f6e2a4))

## [2.1.1](https://github.com/Phenome-Health/kraken/compare/v2.1.0...v2.1.1) (2026-08-18)


### Documentation

* add a contributing.md, other minor tweaks ([16b1568](https://github.com/Phenome-Health/kraken/commit/16b1568201eda07bafbb0d82adcb148cdfca0a7b))
* add zenodo link ([5a6a6ce](https://github.com/Phenome-Health/kraken/commit/5a6a6ce65d9a523d4005d80c6ca6cf07da1eb49a))
* add/adjust badges showing in readme ([7bd7366](https://github.com/Phenome-Health/kraken/commit/7bd73660d2220a34ea8a94a6c3b8efb32787990b))

## [2.1.0](https://github.com/Phenome-Health/kraken/compare/v0.1.0...v2.1.0) (2026-08-15)


### Features

* add bio age ingestion (v1) ([a370f1a](https://github.com/Phenome-Health/kraken/commit/a370f1a6858e77049d7d41ecb7ca09afbe953152))
* add exclude_sources to build config ([66084f5](https://github.com/Phenome-Health/kraken/commit/66084f59df33828c164ca401bb87807a59927496))
* add graph shuffle post-processing step ([29320cc](https://github.com/Phenome-Health/kraken/commit/29320cc2154a902ef2fb2839ef4d83d43e3d4656))
* add initial PGS catalog ingestion (no gene associations yet) ([c2e42af](https://github.com/Phenome-Health/kraken/commit/c2e42afde363f1a48fe430203d2e6ea5e6f17f9f))
* add interactive confirmation of source versions (but not in CI) ([70998c6](https://github.com/Phenome-Health/kraken/commit/70998c6994a6e20df957eafef3a60286342ccc04))
* add joint KLAT to metagraph ([69ce407](https://github.com/Phenome-Health/kraken/commit/69ce407a8bb086d9f4c483db674df8030cc94b9d))
* add KLAT bubble heatmap plot ([cc316c3](https://github.com/Phenome-Health/kraken/commit/cc316c3f3e9f6d7f283f0acd11e50de3e809e8dd))
* add loinc direct ingest ([02841eb](https://github.com/Phenome-Health/kraken/commit/02841ebcc64da33577bbb7787971b10cdcfe97c7))
* add microbiome-kg (and default KS handling) ([7235cce](https://github.com/Phenome-Health/kraken/commit/7235cce1f6cc0074c2c6657b16b8f14697e92c00))
* add multiomics-kg ([c879307](https://github.com/Phenome-Health/kraken/commit/c879307c688cb2f515bc6f0adb6e99be5056a009))
* add NIH CDE ingest ([7df84c6](https://github.com/Phenome-Health/kraken/commit/7df84c6159abb252d9f99d5a1abe1bbd8a5956e3))
* add predicate override mechanism to base harmonizer ([923d05c](https://github.com/Phenome-Health/kraken/commit/923d05c30fc5c793dff7e0305a390108b58f1e9d))
* add primary KS-based edge exclusion mechanism ([f194908](https://github.com/Phenome-Health/kraken/commit/f194908a5ed343dd24b64e361fcdf84adffcddd1))
* add robokop as a source ([fb17f99](https://github.com/Phenome-Health/kraken/commit/fb17f991ee57e3fe1a2ae584bb77e1eea3ce4067))
* add rudimentary subgraph extraction scripts ([27f0bf7](https://github.com/Phenome-Health/kraken/commit/27f0bf71f5b9f0590065c54a192ef3768cb32aa0))
* add script for chord diagram meta-doubles viz ([2c9258c](https://github.com/Phenome-Health/kraken/commit/2c9258c354c14030fc32b476899969f00be3f242))
* add script for comparing metagraph source overlap ([f84984b](https://github.com/Phenome-Health/kraken/commit/f84984bf0c2d05f702ff1ae4a52d8e348a70a985))
* add script for source network visualization; tweak family colors ([42da638](https://github.com/Phenome-Health/kraken/commit/42da638baf30d4942d9e6bc1bb8a230377d977ce))
* add shell for harmonizer validation ([2fb79b4](https://github.com/Phenome-Health/kraken/commit/2fb79b45a82fbc448f0b4d4c302ccd6ff95d2a2b))
* add source versions, biolink version to metagraph ([b0b203a](https://github.com/Phenome-Health/kraken/commit/b0b203a78b6e9a77e7e6d87daec54d98b4e9da65))
* add taxa node prop, skip negated edges ([4876116](https://github.com/Phenome-Health/kraken/commit/487611650dd729259f65fe73e4b6fdaeda1a7c69))
* add translator kg open harmonizer; add trapi source parsing ([4309292](https://github.com/Phenome-Health/kraken/commit/4309292b9057d7c9f1febaf8b6fbfd90bcbe2d88))
* add v1 bio bmi ingestion, some pgs tweaks ([0e971c0](https://github.com/Phenome-Health/kraken/commit/0e971c0c2343cf732654ebac6db9b57fcb1f723a))
* add validation-only option ([c5ae133](https://github.com/Phenome-Health/kraken/commit/c5ae1330374fa8876580deae658a1e41a0000fbc))
* allow multiple dtypes for node/edge schemas ([45012f5](https://github.com/Phenome-Health/kraken/commit/45012f5965feac711e863747fee923967cd05ea5))
* always filter out INCHI curies, uppercase UNIIs ([858245c](https://github.com/Phenome-Health/kraken/commit/858245cbb2c112a27471e5bfef428d1a10511b48))
* auto-fix node/edge types with repeated prefixes ([e0ebc17](https://github.com/Phenome-Health/kraken/commit/e0ebc175cbce71308ece12e2444a60fbeb8f126c))
* build out preliminary validator ([06f2b62](https://github.com/Phenome-Health/kraken/commit/06f2b629820f268dde0366e95578e8289c8b8f32))
* clean up metagraph creation, remove unused bits ([abd12f2](https://github.com/Phenome-Health/kraken/commit/abd12f296e549ceae2fb34ae4a648f26adec5812))
* continue working through molepro overrides ([1ff96bb](https://github.com/Phenome-Health/kraken/commit/1ff96bbd0b626525f3fc0189c61871e1ad5d911f))
* convert rest of harmonizers to new base class paradigm ([d257dac](https://github.com/Phenome-Health/kraken/commit/d257dace9cf0fe61ffc8efba4ebbafda534cedb4))
* derive source_name from source_infores ([08677e7](https://github.com/Phenome-Health/kraken/commit/08677e7382eeb162b6fe86067e9adc1e66a2f5eb))
* exclude ubergraph from robokop ([72e9b55](https://github.com/Phenome-Health/kraken/commit/72e9b554c0dfb05e1db8cc1727a7aeae113b4448))
* extract node urls from robokop, minor cleanup ([a7c7917](https://github.com/Phenome-Health/kraken/commit/a7c791757080890f127babbb181b9c63068bb0c5))
* filter non-leaf categories from merged nodes ([0d11fb7](https://github.com/Phenome-Health/kraken/commit/0d11fb7e0fa512780315b804794a9c17c1e01c3a))
* form tarballs and use class for Orchestrator ([1de6cd9](https://github.com/Phenome-Health/kraken/commit/1de6cd9221795bbc953da60839fde753200c95dd))
* handle improper KLAT of "unspecified" ([5134493](https://github.com/Phenome-Health/kraken/commit/5134493860baba4f51adc025ae5cfbc0083c026e))
* include kraken version in tarball name(s) ([cf3599b](https://github.com/Phenome-Health/kraken/commit/cf3599b644909a73bcdb6562340d06b241b39af2))
* include semmeddb edges from kg2 ([87a36fd](https://github.com/Phenome-Health/kraken/commit/87a36fd3d9bfddfcf6205ccb088cf0068d24adb2))
* make metagraph creation an option vs. step [#4](https://github.com/Phenome-Health/kraken/issues/4) ([d84cd65](https://github.com/Phenome-Health/kraken/commit/d84cd654ad0ccd235e18e92b6f4fb54dcbc33aaa))
* make post-zipping input files a build option ([6d2ae71](https://github.com/Phenome-Health/kraken/commit/6d2ae7134503c1e14f590b9b24c77022d54b5775))
* make publications a top-level node property ([a13e1fd](https://github.com/Phenome-Health/kraken/commit/a13e1fdda9c794fc577a21840a14b2348c73294d))
* merge edges across aggregators (if otherwise match) ([14a1b2e](https://github.com/Phenome-Health/kraken/commit/14a1b2ed974813ef0477250f469c28b7d9819fdd))
* merge edges across aggregators; add disk-based edge merging ([97a9312](https://github.com/Phenome-Health/kraken/commit/97a93121fcd2c5efe296cdafb99517bbdc8e7bf3))
* move to topology-preserving shuffling algorithm ([723f311](https://github.com/Phenome-Health/kraken/commit/723f31133efbc2d79644324786a39bc436866ef1))
* nest meta triples/doubles in metagraph ([2a71969](https://github.com/Phenome-Health/kraken/commit/2a71969a42473ae6afaf7288f1ec8148ed8036c2))
* record graph name/version in metagraph, remove viewer ([a5c892d](https://github.com/Phenome-Health/kraken/commit/a5c892d486822b39b1b0df6e755c74af62a7ee85))
* record supporting source counts w/ primary in metagraphs ([f06c5b0](https://github.com/Phenome-Health/kraken/commit/f06c5b040286394a09d497306d5628943b620dad))
* refine pgs catalog ingest more (edge types, license stuff) ([cb36648](https://github.com/Phenome-Health/kraken/commit/cb36648affcdee4f14d7aa13f1fbff210d6d11e6))
* refine pgs catalog ingestion (add gene/variant edges) ([95103be](https://github.com/Phenome-Health/kraken/commit/95103be3462d6550168d572b817cac887cae8a49))
* report cross-family meta-edges and edge instances ([6ff5563](https://github.com/Phenome-Health/kraken/commit/6ff55635c76a4c61cfcf6f1bcc8f913f17dc4546))
* run all curies through biomapper2's Normalizer ([97b6222](https://github.com/Phenome-Health/kraken/commit/97b62226e75b123812197744115dfea96159e74e))
* separate out qualifiers, detect dynamically ([15c8f0a](https://github.com/Phenome-Health/kraken/commit/15c8f0a454a77d898005efd03982c5d4aceb569e))
* start switch to using base harmonizer, simplify ([f7e5cc0](https://github.com/Phenome-Health/kraken/commit/f7e5cc0fdc72e89d89afa7862d9c1e798d192fd3))
* switch to include_sources list for source selection ([d1ab775](https://github.com/Phenome-Health/kraken/commit/d1ab775b5a6cf8076da95b1129b5916263392440))
* switch to pydantic KrakenConfig! ([bee4b36](https://github.com/Phenome-Health/kraken/commit/bee4b36a01cbe90093ce92f74b3a17447d0b6f54))
* unzip/zip input files before/after use [#3](https://github.com/Phenome-Health/kraken/issues/3) ([fa1b2dc](https://github.com/Phenome-Health/kraken/commit/fa1b2dc6fac40ba8c72d3e783636e47e26f382af))
* write typed build_info.json with schema-drift contract ([5af052c](https://github.com/Phenome-Health/kraken/commit/5af052c8d5a3a882020c1e4c52baf61a4509803e))


### Bug Fixes

* add back 'shuffle_graph' top-level function ([9bfb362](https://github.com/Phenome-Health/kraken/commit/9bfb3622a2587e60027240939e35bc8dc40e769f))
* adjust project root determination post-refactor ([e7f50ff](https://github.com/Phenome-Health/kraken/commit/e7f50ff6b6a6944de198fd92a95716441c9e791a))
* always convert exact mass to float ([17b5321](https://github.com/Phenome-Health/kraken/commit/17b5321918f37147e43b1d5b62ba158508749be5))
* avoid INCHI ids, add category override mechanism ([f7c7ed1](https://github.com/Phenome-Health/kraken/commit/f7c7ed188d2a65a47c76d8ece6979f6a6d023183))
* change bBMI to only ingest trait associations in v1 ([39f0a99](https://github.com/Phenome-Health/kraken/commit/39f0a9925a45666813c8719c87187920c30858b7))
* correct how empty fields are detected, refine list parsing ([7ddbd83](https://github.com/Phenome-Health/kraken/commit/7ddbd8354e805a66de9deb1f15300fa7d62cd7d2))
* create integrated debug dir before use ([42f97ce](https://github.com/Phenome-Health/kraken/commit/42f97ce4404412edff6d07d55b77a5c28ed48867))
* don't load translator equiv ids, instead of can_merge = False ([2472b43](https://github.com/Phenome-Health/kraken/commit/2472b43fd85bfa5255f81c7a0b7260cab157f526))
* dont let translator kg merge existing nodes, other cleanup ([eb63d73](https://github.com/Phenome-Health/kraken/commit/eb63d73c4e7eecae2de88b72deeb020430974242))
* drop 'ignore' props from top-level too (not just attributes) ([f406fd9](https://github.com/Phenome-Health/kraken/commit/f406fd9f447b41e006ba9299404d31e6878f1227))
* equiv ids not being cleared between validator runs ([3e6c23a](https://github.com/Phenome-Health/kraken/commit/3e6c23ae40502a1a55ff07147f1ce3c77bcc2a11))
* exclude supporting data sources from edge keys (causes misses) ([cf5dc58](https://github.com/Phenome-Health/kraken/commit/cf5dc58f107c9da8c517444d02a9d0d9f56ada7a))
* get rid of old (now redundant) organization of attrs by source ([ac5190b](https://github.com/Phenome-Health/kraken/commit/ac5190babeb4cf51d791ceaaca77116740eae943))
* how scripts get harmonized paths post-config changes ([6258a17](https://github.com/Phenome-Health/kraken/commit/6258a17e7e6a49578ce1c1e5e7bdbd22c12c1601))
* log warning (vs. error) on unknown resource_role ([963cb7e](https://github.com/Phenome-Health/kraken/commit/963cb7e06af5d21e1c3133b559faa6e4054c2f55))
* make create node/edge static again, for use in non-base classes ([c41ed32](https://github.com/Phenome-Health/kraken/commit/c41ed3224f53e15ec727b23b0d738ef80e8eb1e5))
* minor bio bmi name/synonym tweaks ([c2a99cc](https://github.com/Phenome-Health/kraken/commit/c2a99cce496aceed4ed0e2558a738d8e21c58937))
* only add attributes dict if not empty ([27ea08d](https://github.com/Phenome-Health/kraken/commit/27ea08d453e79ae3fb497f44c5e17203e4590d8d))
* only pass a few kinds of curies to biomapper (pending bug fixes) ([80106f7](https://github.com/Phenome-Health/kraken/commit/80106f71ba370690fa3290612961344d87e6e6d5))
* pass graph version into test metagraph creation ([4dafe4d](https://github.com/Phenome-Health/kraken/commit/4dafe4df729ccb774cebb74bacbb3ab3f3eeb7dd))
* prevent empty string values in harmonized jsonl, minor fixes ([5eb4fd5](https://github.com/Phenome-Health/kraken/commit/5eb4fd5a72a6e81fe17ed74002cd91bf373e4f11))
* record which fields may be lists, other minor fixes ([69b4a32](https://github.com/Phenome-Health/kraken/commit/69b4a32c64b70c956104a4657a43db31a4f912e3))
* remove publications from HMDB edges (unreliable) ([3672e39](https://github.com/Phenome-Health/kraken/commit/3672e39c0c4f65a4ce5c047625eeadec6b54468b))
* shuffled graph should be written to specified dir ([31ce5af](https://github.com/Phenome-Health/kraken/commit/31ce5af30e8adccd58256bfb8605b1b9a915a7f3))
* specify precise umls version used ([3df2ea9](https://github.com/Phenome-Health/kraken/commit/3df2ea90a8320cf0567246cd17b6753611609521))
* start merged edge count from 0 ([97c9bb3](https://github.com/Phenome-Health/kraken/commit/97c9bb35532966d07c271f08015f1f1fd65ea582))
* tally node prefixes from equiv IDs (not just ID) ([d10a33b](https://github.com/Phenome-Health/kraken/commit/d10a33b655b67f353e765c59487afc0816f23e66))
* temporarily add back aggregator KS to edge keys ([815c007](https://github.com/Phenome-Health/kraken/commit/815c007e4cb0601e31517a83573cb617efe7cc42))
* tiny legend wording fix ([7981989](https://github.com/Phenome-Health/kraken/commit/7981989d17c81065f4c34a91eefe1d227e9d8e92))
* update category viz for kraken 2.0 ([d3701f4](https://github.com/Phenome-Health/kraken/commit/d3701f48907617a222832e87ef283962bb651b88))
* update normalizer usages for new biomapper2 ([4f7ca8e](https://github.com/Phenome-Health/kraken/commit/4f7ca8e061fa4f72c82d3e45131b8440b32b1cce))
* update one more biomapper2 usage ([84cfdff](https://github.com/Phenome-Health/kraken/commit/84cfdfffa7e4dd9c28f8d0b049bad9aa9444f9f7))
* update other external uses of static create_node/edge ([70f7835](https://github.com/Phenome-Health/kraken/commit/70f78358019994ec52e89fe390245208fe83e564))
* use correct config slot names for zipping/unzipping ([3d42c6d](https://github.com/Phenome-Health/kraken/commit/3d42c6d4fffbeb8a839c394589f033fbf27b94d0))
* use valid robokop/refmet infores curies ([043ed88](https://github.com/Phenome-Health/kraken/commit/043ed881cc9f40992e02a4e9afcb754e6c7f0bbe))


### Performance Improvements

* add property name constants for quick access ([d33b3cc](https://github.com/Phenome-Health/kraken/commit/d33b3cc7cb9a79c38b9c07b1a9cf14542b66a93c))
* cache normalized curies for use in edge harmonization ([c717be1](https://github.com/Phenome-Health/kraken/commit/c717be1e022037dc2cd6f6431537647bd194718c))


### Miscellaneous Chores

* release 2.1.0 ([89a4456](https://github.com/Phenome-Health/kraken/commit/89a4456211d20deb297d35e613edcd62339916ec))
