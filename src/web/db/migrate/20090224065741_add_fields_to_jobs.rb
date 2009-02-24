class AddFieldsToJobs < ActiveRecord::Migration
  def self.up
    add_column :jobs, :locus_file_id, :integer
    add_column :jobs, :genotype_file_id, :integer
    add_column :jobs, :results, :text
  end

  def self.down
    remove_column :jobs, :locus_file_id
    remove_column :jobs, :genotype_file_id
    remove_column :jobs, :results
  end
end
