## This code is the source code implementation for the paper "DP-STDP Synthetic Trajectory Data Publishing".



## Abstract

![](/pic/arc.png)

With the popularization of mobile devices with rich sensing capabilities, mobile crowdsensing (MCS) is an emerging paradigm to facilitate urban sensing. In this paper, a synthetic trajectory data release framework DP-STDP based on differential privacy is proposed to protect the user’s trajectory privacy in the scenario where mobile users publish trajectory data to the sensing platform for a period of time. This framework protects the original trajectory by generating synthetic trajectory privacy. First, the geographic space is discretized using location-diversity adaptive grid partition, then the distribution of trajectory spatial characteristics is established to ensure the validity of the data, and an appropriate amount of differentially private noise is added to each part to finally generate synthetic trajectories. At the same time, the impact of Bayesian adversary inference attack is also considered, and trajectory data with high availability is released while protecting the privacy of trajectory data. Finally, experiments show that our method not only provides strong privacy guarantees but also outperforms existing methods in terms of data availability.



## Experimental Environment

**Operating environment：**

AMD Ryzen 7 5800H CPU and 16GB Memory,Ubuntu



**Datasets:** The experiments were conducted using three different datasets: Gowalla, Geolife GPS Trajectories, and a synthetic dataset generated using the Brinkhoff generator.

- **Gowalla:** User trajectory data from the location-based social networking site Gowalla, collected between February 2009 and October 2010, focused on New York City.
- **Geolife GPS Trajectories:** Trajectories collected in Beijing, China, with preprocessing to remove outliers.
- **Brinkhoff:** A synthetic dataset generated for the city of Oldenburg, Germany, using the Brinkhoff generator.



**Evaluation Metrics:**

- **Query Error:** Measures the accuracy of answer count queries on the synthetic dataset relative to the original dataset.
- **Frequent Pattern Average Relative Error:** Assesses the preservation of frequent patterns in the synthetic data.
- **Trip Error:** Evaluates the correlation between the start and end points of trajectories.
- **Diameter Error:** Measures the correlation between the diameters of the original and synthetic trajectories.



**Experimental Settings:**

- The evaluation compared DP-STDP with existing methods (DPT and DP-Star) by varying the size of the differential privacy budget ε (0.1, 0.5, 1.0, and 2.0).
- Each experiment was repeated 10 times, and the results were averaged to ensure reliability.



### Experimental Results

The experimental results in this study demonstrate the effectiveness of the DP-STDP framework in maintaining the privacy and utility of trajectory data in mobile crowdsensing scenarios. The experiments were conducted using three datasets: Gowalla, Geolife GPS Trajectories, and a synthetic dataset generated using the Brinkhoff generator. The evaluation metrics included query error, frequent pattern average relative error, trip error, and diameter error. The results show that the DP-STDP method significantly outperforms existing methods, such as DPT and DP-Star, in terms of data utility while providing strong privacy guarantees. Specifically, the DP-STDP method achieves lower query error and frequent pattern average relative error, indicating better preservation of spatial distribution and movement patterns. Additionally, the trip error and diameter error results confirm the effectiveness of the adaptive grid partition and Markov model in capturing and maintaining the user's mobility patterns and trajectory characteristics.

![](/pic/1.png)

![](/pic/2.png) 
